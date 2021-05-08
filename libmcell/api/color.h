/******************************************************************************
 *
 * Copyright (C) 2021 by
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

#ifndef API_COLOR_H
#define API_COLOR_H

#include "generated/gen_color.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Color: public GenColor {
public:
  COLOR_CTOR()

  void postprocess_in_ctor() override;

  void set_red(const double new_red_) override {
    red = new_red_;
    components_to_rgba();
  }

  void set_green(const double new_green_) override {
    green = new_green_;
    components_to_rgba();
  }

  void set_blue(const double new_blue_) override {
    blue = new_blue_;
    components_to_rgba();
  }

  void set_alpha(const double new_alpha_) override {
    alpha = new_alpha_;
    components_to_rgba();
  }

  void set_rgba(const uint new_rgba_) override {
    rgba = new_rgba_;
    rgba_to_components();
  }
private:
  void components_to_rgba();
  void rgba_to_components();
  void check_component_range(const double value, const char* name);
};

} // namespace API
} // namespace MCell

#endif // API_COLOR_H
