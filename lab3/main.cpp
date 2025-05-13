#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <windows.h>

std::mutex mtx;
std::condition_variable cv;

void worker() {

    std::cout << "Worker: received signal!" << std::endl;
}

void notifier() {

    while (true)
    {
        SHORT keystate = GetKeyState(VK_SPACE);
        if (keystate & 0x8000)
        {
            worker();
        }
        Sleep(30);
    }
}

int main() {

    std::thread t1(notifier);


    //тут уже свободно

    t1.join();
    return 0;
}
