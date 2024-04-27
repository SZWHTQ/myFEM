#pragma once
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

class ThreadPool {
   public:
    explicit ThreadPool(size_t numThreads);

    ~ThreadPool();

    template <class F>
    void enqueue(F&& task) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.emplace(std::forward<F>(task));
        }

        condition.notify_one();
    }

    bool stop;

   private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
};