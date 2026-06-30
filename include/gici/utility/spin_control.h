/**
* @Function: Tick control for spinning
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#include <atomic>
#include <iostream>
#include <vikit/timer.h>
#include <thread>

namespace gici {

// Spin control
class SpinControl {
public:
  SpinControl(double duration);
  ~SpinControl() { }

  // Sleep to ensure spinning rate and restart timer
  void sleep();

  // (Re)set loop duration
  void setDuration(double duration) { duration_ = duration; }

  // Check global status
  static bool ok() { 
    while (wait_.load()) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(int(1e5)));
    }
    return ok_.load();
  }

  // Shutdown all spin controllers
  static void shutdown() {
    ok_.store(false);
    wait_.store(false);
    shutdown_requested_.store(true);
  }

  // Backward-compatible alias for shutdown
  static void kill() { shutdown(); }

  // Request a graceful shutdown that lets the owner stop threads explicitly.
  static void requestShutdown() { shutdown_requested_.store(true); }

  // Check whether a graceful shutdown has been requested.
  static bool shutdownRequested() { return shutdown_requested_.load(); }

  // All spin should wait
  static void wait() { wait_.store(true); }

  // All spin end waiting
  static void run() { wait_.store(false); }

  // Reset controller state for a fresh run
  static void reset() {
    ok_.store(true);
    wait_.store(true);
    shutdown_requested_.store(false);
  }

private:
  // Tick controllers
  vk::Timer timer_;
  double duration_;

  // Whether enable spin in all objects
  static std::atomic<bool> ok_;

  // If all spin should wait
  static std::atomic<bool> wait_;

  // Whether an orderly shutdown has been requested.
  static std::atomic<bool> shutdown_requested_;
};

}
