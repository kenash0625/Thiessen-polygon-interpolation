#include <vector>
#include <string>
#include <map>
#include <stack>
#include <list>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <future>
using namespace std;

struct ac_threads;
enum ac_thread_stat {STOP,RUN};
struct ac_thread
{
	thread t;
	ac_thread_stat s;
	void thread_fun(ac_threads *ths);	
	ac_thread(ac_threads *p);
private:
	ac_thread(const ac_thread &r);
	ac_thread &operator=(const ac_thread &r);
};
struct ac_thread_stat_
{
	ac_thread_stat &m_s;
	ac_thread_stat_(ac_thread_stat &s);
	~ac_thread_stat_();
};
//sender receiver
//使用的人不用知道receiver 只是send
struct ac_threads {
private:
	friend struct ac_thread;
	void bigger(size_t cnt);

	// need to keep track of threads so we can join them
	std::list< ac_thread* > m_workers,m_idles;
	// the task queue
	std::list< std::packaged_task<void()> > m_tasks;

	// synchronization
	std::mutex m_queuemutex,m_thrmutex;
	std::condition_variable m_cv;
	bool m_stop;
public:
	ac_threads(size_t threads);
	template<class F, class... Args>
	auto push(F&& f, Args&&... args)
		//->std::future<typename std::result_of<F(Args...)>::type>
	{
		//using return_type = typename std::result_of<F(Args...)>::type;

		auto task = std::packaged_task<void()>(
			std::bind(std::forward<F>(f), std::forward<Args>(args)...)
			);

		std::future<void> res = task.get_future();
		{
			std::unique_lock<std::mutex> lock(m_queuemutex);
			// don't allow enqueueing after stopping the pool
			if (m_stop)
			{
				cout << "##$$" << endl;
				throw std::runtime_error("enqueue on stopped ThreadPool");
			}
			m_tasks.push_back(std::move(task));// .emplace_back([task]() { (task)(); });
		}
		{
			std::lock_guard<std::mutex> lkg(m_thrmutex);
			auto iter = find_if(m_workers.begin(), m_workers.end(), [](ac_thread *w)->bool {return w->s == STOP; });
			if (iter == m_workers.end())
			{
				bigger(4);
			}
		}
		m_cv.notify_one();
		return res;
	}
	void wait();
	~ac_threads();

};

struct ac_node
{
	string cd;
	string to;
	vector<string> from;
};
struct ac_calibration
{
	//ac_notifylist<ac_msg> resnotify;
	map<string, ac_node> allnode;

	void trace(string &root,string &rt);
	void wait();
	void stringsplit(const string &r, vector<string> &v);
};
