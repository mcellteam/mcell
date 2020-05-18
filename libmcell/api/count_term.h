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

#ifndef API_COUNT_TERM_H
#define API_COUNT_TERM_H

#include "../generated/gen_count_term.h"
#include "../api/common.h"

namespace MCell {
namespace API {

class CountTerm: public GenCountTerm, public std::enable_shared_from_this<CountTerm> {
public:
  COUNT_TERM_CTOR()

  std::shared_ptr<CountTerm> create_expr_term(ExprNodeType op, std::shared_ptr<CountTerm> op2) {
    std::shared_ptr<CountTerm> res = std::make_shared<CountTerm>();
    res->node_type = op;
    res->left_node = shared_from_this();
    res->right_node = op2;
    return res;
  }

  std::shared_ptr<CountTerm> __add__(std::shared_ptr<CountTerm> op2) override {
    return create_expr_term(ExprNodeType::Add, op2);
  }

  std::shared_ptr<CountTerm> __sub__(std::shared_ptr<CountTerm> op2) override {
    return create_expr_term(ExprNodeType::Add, op2);
  }

};

} // namespace API
} // namespace MCell

#endif // API_COUNT_TERM_H
