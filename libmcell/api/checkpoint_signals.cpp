/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include "checkpoint_signals.h"

#include <iostream>
#include <set>
#include <signal.h>

#include "api/common.h"
#include "api/model.h"

using namespace std;

namespace MCell {
namespace API {

// WARNING: not multithread-safe
std::set<Model*> g_models;
struct sigaction g_previous_sigaction_sigusr1;
struct sigaction g_previous_sigaction_sigusr2;
struct sigaction g_previous_sigaction_sigalrm;

void checkpoint_signal_handler(int signo) {
  // checkpoint can be requested multiple times even when the
  // scheduling of a checkpoint is running

  bool continue_simulation = false;

  for (Model* m: g_models) {

    // printout for each model instance
#ifndef _WIN32
    if (signo == SIGUSR1) {
      cout << "User signal SIGUSR1 detected, scheduling a checkpoint and continuing simulation.\n";
      continue_simulation = true;
    }
    else if (signo == SIGUSR2) {
      cout << "User signal SIGUSR2 detected, scheduling a checkpoint and terminating simulation afterwards.\n";
      continue_simulation = false;
    }
#endif
    if (signo == SIGALRM) {
      cout << "Signal SIGALRM detected - periodic or time limit elapsed, scheduling a checkpoint ";
      if (m->config.continue_after_sigalrm) {
        cout << "and continuing simulation.\n";
        continue_simulation = true;
      }
      else {
        cout << "and terminating simulation afterwards.\n";
        continue_simulation = false;
      }
    }

    m->schedule_checkpoint(0, continue_simulation);
  }
}


// Set signal handlers for checkpointing on SIGUSR signals.
void set_checkpoint_signals(Model* model) {

  bool already_set = !g_models.empty();

  g_models.insert(model);

  // Windows does not support USR signals
  if (!already_set) {
#ifndef _WIN32
    struct sigaction sa;
    sa.sa_sigaction = NULL;
    sa.sa_handler = &checkpoint_signal_handler;
    sa.sa_flags = SA_RESTART;
    sigfillset(&sa.sa_mask);

    if (sigaction(SIGUSR1, &sa, &g_previous_sigaction_sigusr1) != 0) {
      throw RuntimeError("Failed to install SIGUSR1 signal handler.");
    }
    if (sigaction(SIGUSR2, &sa, &g_previous_sigaction_sigusr2) != 0) {
      throw RuntimeError("Failed to install SIGUSR2 signal handler.");
    }
#endif
    if (sigaction(SIGALRM, &sa, &g_previous_sigaction_sigalrm) != 0) {
      throw RuntimeError("Failed to install SIGUSR2 signal handler.");
    }
  }
}


void unset_checkpoint_signals(Model* model) {
  if (g_models.count(model) == 0) {
    // either not set or unset twice, ignore
    return;
  }

  g_models.erase(model);

  // Windows does not support USR signals
  if (g_models.empty()) {
#ifndef _WIN32
    if (sigaction(SIGUSR1, &g_previous_sigaction_sigusr1, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGUSR1 signal handler.\n";
    }
    if (sigaction(SIGUSR2, &g_previous_sigaction_sigusr2, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGUSR2 signal handler.\n";
    }
#endif
    if (sigaction(SIGALRM, &g_previous_sigaction_sigalrm, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGALRM signal handler.\n";
    }
  }
}

} // namespace API
} // namespace MCell

