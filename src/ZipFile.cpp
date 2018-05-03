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

#include "ZipFile.h"

#include <zip.h>

#include <errno.h>
#include <stdexcept>


ZipFile::ZipFile(const std::string& filename) {
	if (filename.empty())
		throw std::runtime_error("ZipFile::ZipFile(): Empty file name given");

	int err;
	char errstr[1024];

	// Open ZIP file
	m_zip = zip_open(filename.c_str(), 0, &err);
	if (!m_zip) {
		zip_error_to_str(errstr, sizeof(errstr), err, errno);
		
		std::string msg("ZipFile::ZipFile(): Cannot open zip file ");
		msg.append(filename);
		msg.append(" : ");
		msg.append(errstr);
		
		throw std::runtime_error(msg);
	}

	// Get number of files
	m_nEntries = zip_get_num_files(m_zip); // Mac OS X 10.6
	//m_nEntries = zip_get_num_files(m_zip, 0); // All other OSs
}

ZipFile::~ZipFile() {
	if (m_zip)
		zip_close(m_zip);
}

void ZipFile::loadFile(zip_int64_t index, char*& buffer, zip_int64_t& bufferLength) {
	if (!m_zip)
		throw std::runtime_error("ZipFile::loadFile(): No ZIP file open");

	if (index >= m_nEntries)
		throw std::runtime_error("ZipFile::loadFile(): Invalid file index");

	struct zip_stat stats;

	// Get file size
	if (zip_stat_index(m_zip, index, 0, &stats) != 0) {
		std::string msg("ZipFile::loadFile(): Cannot find file in ZIP: ");
		msg.append(zip_strerror(m_zip));

		throw std::runtime_error(msg);
	}

	bufferLength = stats.size;

	// Open file and allocate memory
	zip_file* f = zip_fopen_index(m_zip, index, 0);
	if (f == NULL) {
		std::string msg("ZipFile::loadFile(): Cannot open file in ZIP: ");
		msg.append(zip_strerror(m_zip));

		throw std::runtime_error(msg);
	}

	buffer = new char[bufferLength];

	// Uncompress file into memory
	if (zip_fread(f, buffer, bufferLength) != bufferLength) {
		std::string msg("ZipFile::loadFile(): Cannot uncompress file in ZIP: ");
		msg.append(zip_file_strerror(f));

		zip_fclose(f);
		throw std::runtime_error(msg);
	}

	zip_fclose(f);
}

