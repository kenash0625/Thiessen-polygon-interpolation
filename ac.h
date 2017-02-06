#include <vector>
#include <string>
#include <map>
#include <list>
#include <mutex>
#include <condition_variable>
using namespace std;
struct ac_msg
{
	std::thread::id tid;
	string cd,rt;
};
template<class T>
class ac_notifylist
{
	std::condition_variable cv;
	std::mutex m;
	std::list<T> stk;
public:
	void pop(T &i)
	{
		std::unique_lock<std::mutex> lk(m);
		cv.wait(lk, [this] {return !stk.empty(); });//满足条件-获取锁 不满足-释放锁 等待被唤醒
		i = *stk.begin();
		stk.pop_front();
	}
	void push(T &i)
	{
		std::lock_guard<std::mutex> lk(m);
		stk.push_back(i);
		cv.notify_one();
	}
};

struct ac_node
{
	string cd;
	string to;
	vector<string> from;
};
struct ac_calibration
{
	ac_notifylist<ac_msg> resnotify;
	map<string, ac_node> allnode;

	int trace(const string &root,const string &rt);
	void wait();
	void stringsplit(const string &r, vector<string> &v);
};
