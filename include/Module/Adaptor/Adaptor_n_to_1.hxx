#include <string>
#include <sstream>

#include "Tools/Exception/exception.hpp"
#include "Module/Adaptor/Adaptor_n_to_1.hpp"

namespace aff3ct
{
namespace module
{

Task& Adaptor_n_to_1
::operator[](const adp::tsk t)
{
	return Module::operator[]((size_t)t);
}

Socket& Adaptor_n_to_1
::operator[](const adp::sck::push_n s)
{
	return Module::operator[]((size_t)adp::tsk::push_n)[(size_t)s];
}

Socket& Adaptor_n_to_1
::operator[](const adp::sck::pull_1 s)
{
	return Module::operator[]((size_t)adp::tsk::pull_1)[(size_t)s];
}

Adaptor_n_to_1
::Adaptor_n_to_1(const size_t n_elmts,
                 const std::type_index datatype,
                 const size_t buffer_size,
                 const bool active_waiting,
                 const int n_frames)
: Adaptor(n_elmts, datatype, buffer_size, n_frames),
  active_waiting(active_waiting),
  cnd_put (new std::vector<std::condition_variable>(1000)),
  mtx_put (new std::vector<std::mutex             >(1000)),
  cnd_pull(new             std::condition_variable (    )),
  mtx_pull(new             std::mutex              (    ))
{
	this->init();
}

Adaptor_n_to_1
::Adaptor_n_to_1(const std::vector<size_t> &n_elmts,
                 const std::vector<std::type_index> &datatype,
                 const size_t buffer_size,
                 const bool active_waiting,
                 const int n_frames)
: Adaptor(n_elmts, datatype, buffer_size, n_frames),
  active_waiting(active_waiting),
  cnd_put (new std::vector<std::condition_variable>(1000)),
  mtx_put (new std::vector<std::mutex             >(1000)),
  cnd_pull(new             std::condition_variable (    )),
  mtx_pull(new             std::mutex              (    ))
{
	this->init();
}

void Adaptor_n_to_1
::init()
{
	const std::string name = "Adaptor_n_to_1";
	this->set_name(name);
	this->set_short_name(name);

	std::function<size_t(Task&, const size_t, const std::type_index&, const std::string&)> create_socket_in =
		[this](Task& p, const size_t n_elmts, const std::type_index& datatype, const std::string& n)
		{
			     if (datatype == typeid(int8_t )) return this->template create_socket_in<int8_t >(p, n, n_elmts);
			else if (datatype == typeid(int16_t)) return this->template create_socket_in<int16_t>(p, n, n_elmts);
			else if (datatype == typeid(int32_t)) return this->template create_socket_in<int32_t>(p, n, n_elmts);
			else if (datatype == typeid(int64_t)) return this->template create_socket_in<int64_t>(p, n, n_elmts);
			else if (datatype == typeid(float  )) return this->template create_socket_in<float  >(p, n, n_elmts);
			else if (datatype == typeid(double )) return this->template create_socket_in<double >(p, n, n_elmts);
			else
			{
				std::stringstream message;
				message << "This should never happen.";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}
		};

	std::function<size_t(Task&, const size_t, const std::type_index&, const std::string&)> create_socket_out =
		[this](Task& p, const size_t n_elmts, const std::type_index& datatype, const std::string& n)
		{
			     if (datatype == typeid(int8_t )) return this->template create_socket_out<int8_t >(p, n, n_elmts);
			else if (datatype == typeid(int16_t)) return this->template create_socket_out<int16_t>(p, n, n_elmts);
			else if (datatype == typeid(int32_t)) return this->template create_socket_out<int32_t>(p, n, n_elmts);
			else if (datatype == typeid(int64_t)) return this->template create_socket_out<int64_t>(p, n, n_elmts);
			else if (datatype == typeid(float  )) return this->template create_socket_out<float  >(p, n, n_elmts);
			else if (datatype == typeid(double )) return this->template create_socket_out<double >(p, n, n_elmts);
			else
			{
				std::stringstream message;
				message << "This should never happen.";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}
		};

	auto &p1 = this->create_task("push_n", (int)adp::tsk::push_n);
	std::vector<size_t> p1s_in;
	for (size_t s = 0; s < this->n_sockets; s++)
		p1s_in.push_back(create_socket_in(p1, this->n_elmts[s], this->datatype[s], "in" + std::to_string(s)));

	this->create_codelet(p1, [p1s_in](Module &m, Task &t) -> int
	{
		std::vector<const int8_t*> sockets_dataptr(p1s_in.size());
		for (size_t s = 0; s < p1s_in.size(); s++)
			sockets_dataptr[s] = static_cast<const int8_t*>(t[p1s_in[s]].get_dataptr());

		static_cast<Adaptor_n_to_1&>(m).push_n(sockets_dataptr);
		return 0;
	});

	auto &p2 = this->create_task("pull_1", (int)adp::tsk::pull_1);
	std::vector<size_t> p2s_out;
	for (size_t s = 0; s < this->n_sockets; s++)
		p2s_out.push_back(create_socket_out(p2, this->n_elmts[s], this->datatype[s], "out" + std::to_string(s)));

	this->create_codelet(p2, [p2s_out](Module &m, Task &t) -> int
	{
		std::vector<int8_t*> sockets_dataptr(p2s_out.size());
		for (size_t s = 0; s < p2s_out.size(); s++)
			sockets_dataptr[s] = static_cast<int8_t*>(t[p2s_out[s]].get_dataptr());

		static_cast<Adaptor_n_to_1&>(m).pull_1(sockets_dataptr);
		return 0;
	});
}
}
}