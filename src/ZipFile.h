/* 
 * CoroEval
 *
 * An evaluation tool for coronary artery reconstructions.
 *
 * Copyright © 2014:
 *
 * Christoph Forman, Universität Erlangen-Nürnberg
 * christoph.forman@cs.fau.de
 *
 * Chris Schwemmer, Universität Erlangen-Nürnberg
 * chris.schwemmer@cs.fau.de
 *
 * Jens Wetzl, Universität Erlangen-Nürnberg
 * jens.wetzl@cs.fau.de
 *
 * CoroEval is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoroEval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __ZIPFILE_H__
#define __ZIPFILE_H__

#include <zipconf.h>

#include <string>


// Forward declarations
struct zip;

// C++ wrapper class for libzip functionality
class ZipFile {
public:
	ZipFile(const std::string& filename);
	~ZipFile();

	zip_int64_t nEntries() const { return m_nEntries; }

  void loadFile(zip_int64_t index, char*& buffer, zip_int64_t &bufferLength);
	
private:
	zip* m_zip;
	zip_int64_t m_nEntries;
};

#endif

