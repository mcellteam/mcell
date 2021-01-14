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
#include <unistd.h>

#include "api/common.h"
#include "api/model.h"
#include "api/python_exporter.h"
#include "src4/world.h"
#include "src4/viz_output_event.h"

using namespace std;

#ifndef _WIN64
struct sigaction g_previous_sigaction_sigusr1;
struct sigaction g_previous_sigaction_sigusr2;
struct sigaction g_previous_sigaction_sigalrm;
#endif

namespace MCell {
namespace API {

// WARNING: not multithread-safe
std::set<Model*> g_models;

// WARNING: only limited set of calls is allowed in signal handlers,
// e.g. no malloc
void checkpoint_signal_handler(int signo) {
  // checkpoint can be requested multiple times even when the
  // scheduling of a checkpoint is running

  for (Model* m: g_models) {
    if (m->get_world() != nullptr) {
      m->get_world()->set_to_create_checkpoint_event_from_signal_hadler(signo, m);
    }
  }
}


// Set signal handlers for checkpointing on SIGUSR signals.
void set_checkpoint_signals(Model* model) {

  bool already_set = !g_models.empty();

  g_models.insert(model);

  // Windows does not support USR signals
  if (!already_set) {
#ifndef _WIN64
    struct sigaction sa;
    sa.sa_sigaction = NULL;
    sa.sa_handler = &checkpoint_signal_handler;
    sa.sa_flags = SA_RESTART;
    sigfillset(&sa.sa_mask);

    if (::sigaction(SIGUSR1, &sa, &g_previous_sigaction_sigusr1) != 0) {
      throw RuntimeError("Failed to install SIGUSR1 signal handler.");
    }
    if (::sigaction(SIGUSR2, &sa, &g_previous_sigaction_sigusr2) != 0) {
      throw RuntimeError("Failed to install SIGUSR2 signal handler.");
    }
    if (::sigaction(SIGALRM, &sa, &g_previous_sigaction_sigalrm) != 0) {
      throw RuntimeError("Failed to install SIGUSR2 signal handler.");
    }
#endif
  }
}


void unset_checkpoint_signals(Model* model) {
  if (g_models.count(model) == 0) {
    // either not set or unset twice, ignore
    return;
  }

  g_models.erase(model);

  // Windows does not support USR signals
  // SIGALRM should be supported somehow but it does not work yet
  if (g_models.empty()) {
#ifndef _WIN32
    if (sigaction(SIGUSR1, &g_previous_sigaction_sigusr1, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGUSR1 signal handler.\n";
    }
    if (sigaction(SIGUSR2, &g_previous_sigaction_sigusr2, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGUSR2 signal handler.\n";
    }
    if (sigaction(SIGALRM, &g_previous_sigaction_sigalrm, nullptr) != 0) {
      cout << "Warning: failed to uninstall SIGALRM signal handler.\n";
    }
#endif
  }
}


void save_checkpoint_func(const float_t time, CheckpointSaveEventContext ctx) {

  const World* world = ctx.model->get_world();

  release_assert(
      world->scheduler.get_event_being_executed()->type_index == EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT &&
      "May be called only from event with index EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT, "
      " world/model data may be inconsistent otherwise"
  );

  uint64_t current_it = world->stats.get_current_iteration();

  std::string dir;
  if (ctx.append_it_to_dir) {
    dir =
        ctx.dir_prefix +
        VizOutputEvent::iterations_to_string(world->stats.get_current_iteration(), ctx.model->config.total_iterations) +
        BNG::PATH_SEPARATOR;
  }
  else {
    dir = ctx.dir_prefix;
  }

  cout << "Saving scheduled checkpoint in iteration " << current_it << " into " << dir << "\n";

  PythonExporter exporter(ctx.model);
  exporter.save_checkpoint(dir);
}

} // namespace API
} // namespace MCell

