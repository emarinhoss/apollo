// WarpX lib includes
#include "wxiobase.h"
#include "wxiotmpl.h"

// std includes
#include <sstream>

WxIoBase::WxIoBase() {
    dumpNo = -1;
}

WxIoBase::WxIoBase(const std::string& sfx) {
    dumpNo = -1;
    suffix = sfx;
}

WxIoBase::WxIoBase(const std::string& bn, int d, const std::string& sfx) {
    baseName = bn;
    dumpNo = d;
    suffix = sfx;
}

WxIoBase::~WxIoBase() {
}

WxIoNodeType WxIoBase::createDumpFile(const std::string& dataName) {
    return createFile(getDumpFileName(dataName));
}

WxIoNodeType WxIoBase::openDumpFile(const std::string& dataName,
                                    const std::string& perms) {
    return openFile(getDumpFileName(dataName), perms);
}

std::string WxIoBase::getDumpFileName(const std::string& dataName) {
    std::stringstream sstr;
    if (baseName.size()) sstr << baseName << "_";
    sstr << dataName;
    if (dumpNo >=0) sstr << "_" << dumpNo;
    sstr << suffix;
    std::string fileName;
    sstr >> fileName;
    return fileName;
}

void WxIoBase::addOpenFile(WxIoNodeType node) {
    openFiles.push_back(node);
}

void WxIoBase::removeOpenFile(WxIoNodeType node) {
    std::vector<WxIoNodeType>::iterator itr;
    for (itr = openFiles.begin(); itr != openFiles.end(); ++itr) {
        if (**itr == *node) {
            openFiles.erase(itr);
            break;
        }
    }
}

void WxIoBase::closeOpenFiles() {
    std::vector<WxIoNodeType>::iterator itr;
    for (itr = openFiles.begin(); itr != openFiles.end(); ++itr) {
        delete *itr;
    }
    openFiles.clear();
}

