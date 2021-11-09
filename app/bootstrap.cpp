#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "aff3ct.hpp"
using namespace aff3ct;

struct params
{
    int M = 256;                     // number of SCSs
    int N = 14;                      // number of time domain symbols

    int scs = 15000;                 // Subcarrier spacing
    int nFFT = 256;
    int nSCS = 256;                  // number of SCS
    float cpRatio = 0.07;            // CP ratio for OFDM

    int crcLength = 24;              // CRC length
    int QAM = 16;                    // QAM modulation order

    int numFrames = 10000;           // number of frame errors
    int seed = 0;                    // PRNG seed for the AWGN channel
    float ebn0_min = 0.00f;          // minimum SNR value
    float ebn0_max = 10.00f;         // maximum SNR value
    float ebn0_step = 1.00f;         // SNR step
    int fe = 100;                    // maximum errors count

    float R = 0.5;                   // code rate (R=K/N)

    int sample_rate = 3840000;          // sample rate
    int bitPerSym;                   // bits per symbol
    int blockSize;
    int infoLength;
    
    std::vector<int> cp;             // CP length of OFDM
};

void init_params(params &p)
{
    p.bitPerSym = log2(p.QAM);
    p.blockSize = p.M * p.N * p.bitPerSym;
    p.infoLength = p.blockSize * p.R;

    int s_cp = 0;
    p.cp.resize(p.N);
    for (auto i = 0; i < p.N; i++)
    {
        p.cp[i] = round(p.nFFT * p.cpRatio);
        s_cp += p.cp[i];
    } 

    int points1s = (s_cp+p.nFFT*p.N)*1000; // sampling points in 1ms
    assert(points1s <= p.sample_rate);
    if (points1s != p.sample_rate)
    {
        std::cout<<"Warning: the extra CPs are added to the 1st symbol to align the Sampling rate."<<std::endl;
        p.cp[0] += p.sample_rate/1000-(s_cp+p.nFFT*p.N);
    }
    
    std::cout << "# * Simulation parameters: " << std::endl;
    std::cout << "#    ** Frames in total  = " << p.numFrames << std::endl;
    std::cout << "#    ** Noise seed       = " << p.seed << std::endl;
    std::cout << "#    ** Info. bits (K)   = " << p.infoLength << std::endl;
    std::cout << "#    ** Frame size       = " << p.blockSize << std::endl;
    std::cout << "#    ** Code rate  (R)   = " << p.R << std::endl;
    std::cout << "#    ** EbN0 min   (dB)  = " << p.ebn0_min << std::endl;
    std::cout << "#    ** EbN0 max   (dB)  = " << p.ebn0_max << std::endl;
    std::cout << "#    ** EbN0 step  (dB)  = " << p.ebn0_step << std::endl;
    std::cout << "#" << std::endl;
}

struct modules
{
    std::unique_ptr<module::Source_random<>> source;
    std::unique_ptr<module::Encoder_repetition_sys<>> encoder;
    std::unique_ptr<tools::Constellation_QAM<>> cnstl;
    std::unique_ptr<module::Modem<>> modem;
    std::unique_ptr<module::Otfs<>> waveform;
    std::unique_ptr<module::Channel<>> channel;
    std::unique_ptr<module::Decoder_repetition_std<>> decoder;
    std::unique_ptr<module::Monitor_BFER<>> monitor;
};

void init_modules(const params &p, modules &m)
{
    m.source = std::unique_ptr<module::Source_random<>>(new module::Source_random<>(p.infoLength));
    m.encoder = std::unique_ptr<module::Encoder_repetition_sys<>>(new module::Encoder_repetition_sys<>(round(1/p.R),p.infoLength));
    m.cnstl = std::unique_ptr<tools::Constellation_QAM<>>(new tools::Constellation_QAM<>(p.bitPerSym));
    m.modem = std::unique_ptr<module::Modem_generic<>>(new module::Modem_generic<>(p.QAM,*m.cnstl));
    m.waveform = std::unique_ptr<module::Otfs<>>(new module::Otfs<>(p.nSCS,p.cp,p.N));
    m.channel = std::unique_ptr<module::Channel_AWGN_LLR<>>(new module::Channel_AWGN_LLR<>(p.N));
    m.decoder = std::unique_ptr<module::Decoder_repetition_std<>>(new module::Decoder_repetition_std<>(round(1/p.R),p.infoLength));
    m.monitor = std::unique_ptr<module::Monitor_BFER<>>(new module::Monitor_BFER<>(p.infoLength, p.fe, p.numFrames));
};

struct buffers
{
    std::vector<int> ref_bits;
    std::vector<int> enc_bits;
    std::vector<float> symbols;
    std::vector<float> noisy_symbols;
    std::vector<float> LLRs;
    std::vector<int> dec_bits;
};

void init_buffers(const params &p, buffers &b)
{
    int cnstl_sym = p.blockSize/p.bitPerSym;
    b.ref_bits = std::vector<int>(p.infoLength);
    b.enc_bits = std::vector<int>(p.blockSize);
    b.symbols = std::vector<float>(cnstl_sym);
    b.noisy_symbols = std::vector<float>(cnstl_sym);
    b.LLRs = std::vector<float>(p.blockSize);
    b.dec_bits = std::vector<int>(p.infoLength);
}

struct utils
{
    std::unique_ptr<tools::Sigma<>> noise;                   // a sigma noise type
    std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
    std::unique_ptr<tools::Terminal_std> terminal;           // manage the output text in the terminal
};

void init_utils(const modules &m, utils &u)
{
    // create a sigma noise type
    u.noise = std::unique_ptr<tools::Sigma<>>(new tools::Sigma<>());
    // report the noise values (Es/N0 and Eb/N0)
    u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_noise<>(*u.noise)));
    // report the bit/frame error rates
    u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*m.monitor)));
    // report the simulation throughputs
    u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*m.monitor)));
    // create a terminal that will display the collected data from the reporters
    u.terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(u.reporters));
}

int main(int argc, char **argv)
{
    // get the AFF3CT version
    const std::string v = "v" + std::to_string(tools::version_major()) + "." +
                          std::to_string(tools::version_minor()) + "." +
                          std::to_string(tools::version_release());

    std::cout << "#----------------------------------------------------------" << std::endl;
    std::cout << "# This is a program using the AFF3CT library (" << v << ")" << std::endl;
    std::cout << "# Feel free to improve it as you want to fit your needs." << std::endl;
    std::cout << "#----------------------------------------------------------" << std::endl;
    std::cout << "#" << std::endl;

    params p;
    init_params(p); // create and initialize the parameters defined by the user
    modules m;
    init_modules(p, m); // create and initialize the modules
    buffers b;
    init_buffers(p, b); // create and initialize the buffers required by the modules
    utils u;
    init_utils(m, u); // create and initialize the utils

    // display the legend in the terminal
    u.terminal->legend();

    // loop over the various SNRs
    for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
    {
        // compute the current sigma for the channel noise
        const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
        const auto sigma = tools::esn0_to_sigma(esn0);

        u.noise->set_values(sigma, ebn0, esn0);

        // update the sigma of the modem and the channel
        //m.modem  ->set_noise(*u.noise);
        //m.channel->set_noise(*u.noise);

        // display the performance (BER and FER) in real time (in a separate thread)
        u.terminal->start_temp_report();

        const float chl_p = 0;
        // run the simulation chain
        while (!m.monitor->frame_limit_achieved() && !u.terminal->is_interrupt())
        {
            m.source->generate(b.ref_bits);
            m.encoder->encode(b.ref_bits.data(), b.enc_bits.data());
            m.modem->modulate(b.enc_bits.data(), b.symbols.data());
            //m.channel->add_noise   (b.symbols,       b.noisy_symbols);
            m.channel->add_noise(b.symbols.data(),b.symbols.data(),b.noisy_symbols.data());
            m.modem->demodulate(b.noisy_symbols.data(),b.noisy_symbols.data(), b.LLRs.data());
            m.decoder->decode_siho(b.LLRs, b.dec_bits);
            m.monitor->check_errors(b.dec_bits, b.ref_bits);
        }

        // display the performance (BER and FER) in the terminal
        u.terminal->final_report();

        // reset the monitor for the next SNR
        m.monitor->reset();
        u.terminal->reset();

        // if user pressed Ctrl+c twice, exit the SNRs loop
        if (u.terminal->is_over())
            break;
    }

    std::cout << "# End of the simulation" << std::endl;

    return 0;
}