#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

// Based on code from quadfoam project

namespace FileDialog {

// -----------------------------------------------------------------------------

std::string openFileName(const std::string &defaultPath = "./.*",
	const std::vector<std::string> &filters = {}, const std::string &desc = "");

std::string saveFileName(const std::string &defaultPath = "./.*",
	const std::vector<std::string> &filters = {}, const std::string &desc = "");

// -----------------------------------------------------------------------------

} // namespace FileDialog
