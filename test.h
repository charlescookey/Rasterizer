#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <future>
#include <functional>
#include <atomic>
#include <condition_variable>

class ThreadPool {
public:
    ThreadPool(size_t numThreads);
    ~ThreadPool();

    template <typename Func, typename... Args>
    auto enqueue(Func&& func, Args&&... args) -> std::future<decltype(func(args...))>;

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    std::atomic<bool> stop;
};

ThreadPool::ThreadPool(size_t numThreads) : stop(false) {
    for (size_t i = 0; i < numThreads; ++i) {
        workers.emplace_back([this]() {
            while (true) {
                std::function<void()> task;

                { // Lock scope
                    std::unique_lock<std::mutex> lock(queueMutex);
                    condition.wait(lock, [this]() { return stop || !tasks.empty(); });

                    if (stop && tasks.empty()) return;

                    task = std::move(tasks.front());
                    tasks.pop();
                }

                task(); // Execute task
            }
            });
    }
}

ThreadPool::~ThreadPool() {
    stop = true;
    condition.notify_all();
    for (std::thread& worker : workers)
        worker.join();
}

template <typename Func, typename... Args>
auto ThreadPool::enqueue(Func&& func, Args&&... args) -> std::future<decltype(func(args...))> {
    using returnType = decltype(func(args...));

    auto task = std::make_shared<std::packaged_task<returnType()>>(
        std::bind(std::forward<Func>(func), std::forward<Args>(args)...)
    );

    std::future<returnType> result = task->get_future();
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        tasks.emplace([task]() { (*task)(); });
    }

    condition.notify_one();
    return result;
}