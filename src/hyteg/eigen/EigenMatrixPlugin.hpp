/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

hyteg::PointND< _Scalar, _Rows > mul( const hyteg::PointND< _Scalar, _Cols >& rhs ) const
{
   hyteg::PointND< _Scalar, _Rows > out;
   out.vector_ = this->derived() * rhs.vector_;
   return out;
}

Eigen::Matrix< _Scalar, _Cols, 1 > mul( const Eigen::Matrix< _Scalar, _Cols, 1 >& rhs ) const
{
   return this->derived() * rhs;
}

void setAll( _Scalar scalar)
{
   this->setConstant( scalar );
}