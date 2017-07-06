#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <sstream>

#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif


ostream &operator<<(ostream &os, axis_type axis)
{
    os << to_string(axis);
    return os;
}


string axis_type_to_string(int axis)
{
    if (axis == AXIS_FREQ)
	return "AXIS_FREQ";
    else if (axis == AXIS_TIME)
	return "AXIS_TIME";
    else if (axis == AXIS_NONE)
	return "AXIS_NONE";

    throw runtime_error("rf_pipelines: internal error: bad argument to axis_type_to_string()");
}


bool file_exists(const string &filename)
{
    struct stat s;

    int err = stat(filename.c_str(), &s);
    if (err >= 0)
        return true;
    if (errno == ENOENT)
        return false;

    throw runtime_error(filename + ": " + strerror(errno));
}


vector<string> listdir(const string &dirname)
{
    vector<string> filenames;

    DIR *dir = opendir(dirname.c_str());
    if (!dir)
	throw runtime_error(dirname + ": opendir() failed: " + strerror(errno));

    ssize_t name_max = pathconf(dirname.c_str(), _PC_NAME_MAX);
    name_max = min(name_max, (ssize_t)4096);

    vector<char> buf(sizeof(struct dirent) + name_max + 1);
    struct dirent *entry = reinterpret_cast<struct dirent *> (&buf[0]);
    
    for (;;) {
	struct dirent *result = nullptr;

	int err = readdir_r(dir, entry, &result);	
	if (err)
	    throw runtime_error(dirname + ": readdir_r() failed");
	if (!result)
	    break;

	filenames.push_back(entry->d_name);
    }

    return filenames;
}


// FIXME currently only creates the last directory, e.g. if called with dirname="/a/b/c"
// it assumes that /a/b exists and only creates directory "c".
void makedirs(const string &dirname)
{
    int err = mkdir(dirname.c_str(), 0777);

    if (!err)
	return;

    if (errno != EEXIST) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
    
    struct stat s;
    err = stat(dirname.c_str(), &s);
    
    if (err < 0) {
	stringstream ss;
	ss << "couldn't stat file " << dirname << ": " << strerror(errno);
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }

    if (!S_ISDIR(s.st_mode)) {
	stringstream ss;
	ss << "couldn't create directory " << dirname << ": file already exists and is not a directory";
	
	string err_msg = ss.str();
	cerr << err_msg << "\n";
	throw runtime_error(err_msg);
    }
}


}  // namespace rf_pipelines
