#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>

#include "myIO.h"

using namespace std;
using namespace myiolib;


//class jFileSys
jFileSys::jFileSys()
{      
    init();

    if(data_dir.empty()) 
    {
//        printf("Could not find Database directory within $HOME !! Fatal error\n"); exit(0);
    }
    home_dir += "/";
    curr_dir += "/";
    data_dir += "/";
    for (int i=0; i<data_subd.size();i++) data_subd[i] = data_dir + data_subd[i] + "/";

}

jFileSys::jFileSys(string name)
{      
    casename = name;
    
    init();
    
//    if(data_dir.empty())
//    {
//        printf("Could not find Database directory within $HOME !! Fatal error\n"); exit(0);
//    }
    if(jio.checkDirectory(casename)) 
    {   printf("found directory corresponding to %s\n",casename.c_str());
        case_dir = curr_dir+"/"+casename;
    }
    case_dir += "/";
    home_dir += "/";
    curr_dir += "/";
//    data_dir += "/";
//    for (int i=0; i<data_subd.size();i++) data_subd[i] = data_dir + data_subd[i] + "/";
    
    // binary files
//    thrmfile = case_dir+"Thermo.data";


}

jFileSys::~jFileSys() 
{  
}

void jFileSys::init()
{  
    jio =*(new JIO());    
    
    struct passwd *pw = getpwuid(getuid());
    char* hd = pw->pw_dir;
    home_dir = *(new string(hd));
    
    char* wd = get_current_dir_name();
    curr_dir = *(new string(wd));
    free(wd);

//    DIR *dp;
//    struct dirent *dirp;
//
//    jString jcurd = *(new jString(curr_dir));
//    if(jcurd.contains("Database"))
//    {
//        int sda = jcurd.lastIndexOf("Database");
//        data_dir = curr_dir.substr(0,sda+8);
//    }
//
//    if(data_dir.empty())
//    {
//        int sep = curr_dir.find_last_of("/\\");
//        string upperD = curr_dir.substr(0,sep+1);
//        while(sep > 0 && jString::contains(upperD,home_dir))
//        {
//            if((dp  = opendir(upperD.c_str())) == NULL)
//            {
//                cout << "Error(" << errno << ") opening " << upperD << endl; abort();
//            }
//
//            vector<string> dirs;
//            while ((dirp = readdir(dp)) != NULL)
//            {
//                string upperF = upperD+string(dirp->d_name);
//                dirs.push_back(upperF);
//                if(!jString::contains(upperF,"Database")) continue;
//                // check if it is a directory
//                if(JIO::is_Directory(upperF) )
//                {
//                    data_dir = upperF; break;
//                }
//            }
//            closedir(dp);
//            if(data_dir.empty())
//            {
//                upperD = upperD.substr(0,sep-1);
//                sep = upperD.find_last_of("/\\");
//                upperD = upperD.substr(0,sep+1);
//            }
//            else break;
//        }
//
//    }
//
//    if(!data_dir.empty())
//    {
//        //changing directory (temporarily)
//        chdir(data_dir.c_str());
//        if((dp  = opendir(data_dir.c_str())) == NULL)
//        {
//            cout << "Error(" << errno << ") opening " << data_dir << endl;  abort();
//        }
//        data_subd = *(new vector<string>());
//        while ((dirp = readdir(dp)) != NULL)
//        {
//            string upperF = string(dirp->d_name);
//            if(JIO::is_Directory(upperF) && upperF[0] != '.') data_subd.push_back(upperF);
//        }
//        closedir(dp);
//        //changing directory back
//        chdir(curr_dir.c_str());
//    }

}

void jFileSys::display()
{
    if(case_dir.empty())
    {
        printf(" no case-directory specified\n");
    }
    else
    {
        printf("case-directory = %s\n",case_dir.c_str());
    }
    printf("home-directory = %s\n",home_dir.c_str());
    printf("work-directory = %s\n",curr_dir.c_str());
    printf("data-directory = %s\n",data_dir.c_str());
    printf("data/sub-directories:\n");
    for (int i=0; i<data_subd.size();i++) printf("\t%s\n",data_subd[i].c_str());
}


//class JIO
bool JIO::exists(string filename)
{
    struct stat info;
    int ret = stat(filename.c_str(),&info);
    if(ret==0) return true;
    else       return false;
}

void JIO::open(string filename, string type)
{
    f = fopen(filename.c_str(),type.c_str()); rewind(f);  opened = true;
}
void JIO::open(const char* filename, const char* type)
{
    f = fopen(filename,type); rewind(f); opened = true;
}
void JIO::close() { if(f != NULL)  {  if(opened) { fflush(f); fclose(f); } f = NULL; } opened=false; }


void JIO::writeInt   (int*    pi) { fwrite(pi,sizeof(int)   ,1,f);  }
void JIO::writeFloat (float*  pr) { fwrite(pr,sizeof(float) ,1,f);  }
void JIO::writeDouble(double* pr) { fwrite(pr,sizeof(double),1,f);  }
void JIO::writeBool  (bool*   pb) { fwrite(pb,sizeof(bool)  ,1,f);  }
void JIO::writeString(string s)
{   const char* cs = s.c_str();
    fwrite(cs,sizeof(char),s.length(),f); fputc(0,f);
}
void JIO::writeInt   (int    i) { fwrite(&i,sizeof(int)   ,1,f);  }
void JIO::writeFloat (float  r) { fwrite(&r,sizeof(float) ,1,f);  }
void JIO::writeDouble(double r) { fwrite(&r,sizeof(double),1,f);  }
void JIO::writeBool  (bool   b) { fwrite(&b,sizeof(bool)  ,1,f);  }

int    JIO::readInt()    { int i   ; fread(&i,sizeof(int)   ,1,f); return i; }
float  JIO::readFloat()  { float r ; fread(&r,sizeof(float) ,1,f); return r; }
double JIO::readDouble() { double r; fread(&r,sizeof(double),1,f); return r; }
bool   JIO::readBool()   { bool b  ; fread(&b,sizeof(bool)  ,1,f); return b; }
string JIO::readString()
{
    int i=0;
    int n=0;
    long int pf0 = ftell(f);
    while((i=fgetc(f)) != 0) n++;
    char* cs = new char[n+1];
    fseek(f,pf0,SEEK_SET);
    for(i=0;i<=n;i++) cs[i] = (char)fgetc(f); //cout << "cs[" << i << "]=" << cs[i] << "\n";
    string s(cs);
    delete [] cs;
    return s;
}

void JIO::writeIntArray   (int N, int*    pa) { fwrite(pa,sizeof(int)   ,N,f); }
void JIO::writeFloatArray (int N, float*  pa) { fwrite(pa,sizeof(float) ,N,f); }
void JIO::writeDoubleArray(int N, double* pa) { fwrite(pa,sizeof(double),N,f); }

int*    JIO::readIntArray    (int N){    int* pa= new int  [N];fread(pa,sizeof(int)   ,N,f);return pa;}
float*  JIO::readFloatArray  (int N){  float* pa= new float[N];fread(pa,sizeof(float) ,N,f);return pa;}
double* JIO::readDoubleArray (int N){ double* pa=new double[N];fread(pa,sizeof(double),N,f);return pa;}

void JIO::dump()
{
    while(!feof(f)) cout << fgetc(f) << "\n";
}


int JIO::makeDirectory(string sdir)
{
    int mkd = 0;
    const char* cdir = sdir.c_str();
    struct stat st;
    if(stat(cdir,&st) == 0)
    {
        if(!S_ISDIR(st.st_mode))
        {
            cout << "error: "+sdir+" exists but is not a directory file!\n";
            if(S_ISREG(st.st_mode)) cout << "It is a regular file\n";
            return -1;
        }
    }
    else
    {
        mkd = mkdir(cdir,S_IRWXU);
    }
    return mkd;
}

int JIO::checkDirectory(string sdir)
{
    const char* cdir = sdir.c_str();
    struct stat st;
    if(stat(cdir,&st) == 0)
    {
        if(!S_ISDIR(st.st_mode))
        {
            cerr << "error: "+sdir+" exists but is not a directory file!\n";
            if(S_ISREG(st.st_mode)) cout << "It is a regular file\n";
            exit(0);
        }
    }
    else
    {
        cerr << "error: "+sdir+" does not exist!\n"; exit(0);
    }
    return (1);
}

bool JIO::is_Directory(string sdir)
{
    bool is_dir = false;
    struct stat st;
    if(stat(sdir.c_str(),&st) == 0)
    {
        if(S_ISDIR(st.st_mode)) is_dir = true;
    }
    return (is_dir);
}
    
bool JIO::is_File(string sdir)
{
    bool is_reg = false;
    struct stat st;
    if(stat(sdir.c_str(),&st) == 0)
    {
        if(S_ISREG(st.st_mode)) is_reg = true;
    }
    return (is_reg);
}

// class jString
jString jString::clean()
{   return clean(" ,;:=#$\t\n\r\f");  }

jString jString::clean(const string separ)
{
	jString empty_str;
	if(s.length() <= 0) { return empty_str; }
    string cleaned = s;
    int beg, loc1, loc2 = 0;
    bool comment = false;
    if(cleaned.compare(0,1,"#")==0) comment = true;
    if(cleaned.length() > 1)
    {
        int beg = cleaned.find_first_not_of(' ',0);
        if(cleaned.compare(beg,2,"//")==0) comment = true;
    }
    if(comment) { return empty_str; }
    // remove trailing comments
    beg = cleaned.find_first_of('/',0);
    if(beg > 0) 
    {
        if(cleaned.compare(beg,2,"//")==0) cleaned.erase(beg,cleaned.length());
    }
    beg = cleaned.find_first_of('#',0);
    if(beg > 0) cleaned.erase(beg,cleaned.length());
    // remove leading and trailing separators
    loc1 = cleaned.find_first_not_of(separ);
    loc2 = cleaned.find_last_not_of( separ);
    cleaned = cleaned.substr(loc1,loc2-loc1+1);
    return jString(cleaned);
} //end clean()


vector<jString> jString::parse()
{   return parse(" ,;:=#$\t\n\r\f");  }

vector<jString> jString::parse(const string sep)
{
    string str = (this->clean(sep)).str();
    vector<jString> tokens;
    int nch = str.length();
    if(nch==0) return tokens;

    int beg = 0;
    int end =-1;
    while(beg < nch)
    {
        beg = str.find_first_not_of(sep,beg);
        if(beg != string::npos && beg >= 0)
        {
            end = str.find_first_of(sep,beg+1);
            if(end==string::npos) end = nch;
            tokens.push_back( jString(str.substr(beg,end-beg)));
            beg = end+1;
        }
        else { beg = nch; }
    }
    return tokens;
}

jString jString::toUpperCase()
{
    int N = s.length();
    char* str = (char *)s.data();
    char* tmp = new char[N];
    int i=0;
    while(str[i]) { tmp[i]=toupper(str[i]); i++; }
    string u(tmp,N);
    delete [] tmp;
    return jString(u);
}

jString jString::whiteOut(const string r)
{
    int N = s.length();
    int M = r.length();
    char* str = (char *)s.data();
    char* rem = (char *)r.data();
    char* tmp = new char[N];
    int i=0;
    while(str[i])
    {
        tmp[i] = str[i];
        for(int j=0;j<M;j++) if(str[i]==r[j]) tmp[i]=' ';
        i++;
    }
    string u(tmp,N);
    delete [] tmp;
	return jString(u);
}

char jString::charAt(int n) { return s[n]; }

jString jString::trim()
{   // removes  trailing blanks
    int loc = s.find_last_not_of(' ');
    string u = s.substr(0,loc+1);
	return jString(u);
}

bool jString::contains(string match)
{
    bool ret = false;
    if(s.find(match) < string::npos) ret=  true;
    return ret;
}

bool jString::matches(int pos, string match)
{
    string subset = s.substr(pos,match.length());
    return (subset==match);
}

bool jString::equals(jString c) { return (s==c.str()); }
bool jString::equals( string c) { return (s==c); }
bool jString::startsWith(string c) { return (s.substr(0,c.length())==c); }
bool jString::endsWith(string c) { int beg =  s.length()-c.length(); return (s.substr(beg,s.length())==c); }

jString jString::concat(jString a)  {  return jString(s+a.s); }
jString jString::concat( string a)  {  return jString(s+a  ); }
jString jString::substring(int beg)          { int end = s.length(); return jString(s.substr(beg,end-beg)); }
jString jString::substring(int beg, int end) {                       return jString(s.substr(beg,end-beg)); }

void jString::replace(const char   c, const char   r) {  s = jString::string_replace(s,c,r); }
void jString::replace(const string c, const string r) {  s = jString::string_replace(s,c,r); }

int jString::indexOf(const char c)      {  return s.find_first_of(c); }
int jString::indexOf(const string c)    
{  
    int pos=-1;
    int beg=0;
    int end=beg+c.length();
    while(end < s.length())
    {
        if(s.substr(beg,end-beg) != c) { beg++; end++; }
        else { pos = beg; break; }
    }
    return pos; 
}

int jString::lastIndexOf (const char c) {  return s.find_last_of(c);  }
int jString::lastIndexOf(const string c)    
{  
    int pos=-1;
    int beg= 0;
    int len =s.length()-c.length();
    while(pos < len)
    {
        int end=beg+c.length();
        int rep = -1;
        while(end < s.length())
        {
            if(s.substr(beg,end-beg) != c) { beg++; end++; }
            else { rep = beg; break; }
        }
        if(rep > pos) 
        {
            pos=rep; beg = end;
        }
        else break;
    }
    return pos; 
}

int    jString::iValue() { return atoi(s.c_str()); }
double jString::dValue() { return atof(s.c_str()); }
float  jString::fValue() { return (float)atof(s.c_str()); }

// static methods
string jString::toString(int value)
{
    ostringstream outstr;
    int v = value;
    int n = 0;
    while(v/=10) { n++; }
    outstr << setw(n) << setfill(' ') << value;
    return outstr.str();
}

string jString::toString(int value, int width)
{
    ostringstream outstr;
    outstr << setw(width) << setfill(' ') << value;
    return outstr.str();
}

string jString::toString(int value, int width, char filler)
{
    ostringstream outstr;
    outstr << setw(width) << setfill(filler) << value;
    return outstr.str();
}

string jString::toString(float value)
{
    ostringstream outstr;
    outstr.setf(ios_base::right, ios_base::adjustfield);
    outstr.setf(ios_base::fixed, ios_base::floatfield);
    //outstr.precision(prec);
    outstr << value;
    return outstr.str();
}

string jString::toString(float value, int prec)
{
    ostringstream outstr;
    outstr.setf(ios_base::right, ios_base::adjustfield);
    outstr.setf(ios_base::fixed, ios_base::floatfield);
    outstr.precision(prec);
    outstr << value;
    return outstr.str();
}

string jString::toString(double value)
{
    ostringstream outstr;
    outstr.setf(ios_base::right, ios_base::adjustfield);
    outstr.setf(ios_base::scientific, ios_base::floatfield);
    outstr << value;
    return outstr.str();
}

string jString::toString(double value, int prec)
{
    ostringstream outstr;
    outstr.setf(ios_base::right, ios_base::adjustfield);
    outstr.setf(ios_base::scientific, ios_base::floatfield);
    outstr.precision(prec);
    outstr << value;
    return outstr.str();
}

vector<string> jString::string_parse(const string sin)
{
    return string_parse(sin," ,;:=#$\t\n\r\f");
}
vector<string> jString::string_parse(const string sin, const string sep)
{
    string str = string_clean(sin,sep);
    vector<string> tokens;
    int nch = str.length();
    if(nch==0) return tokens;

    int beg = 0;
    int end =-1;
    while(beg < nch)
    {
        beg = str.find_first_not_of(sep,beg);
        if(beg != string::npos && beg >= 0)
        {
            end = str.find_first_of(sep,beg+1);
            if(end==string::npos) end = nch;
            tokens.push_back( str.substr(beg,end-beg));
            beg = end+1;
        }
        else { beg = nch; }
    }
    return tokens;
}

string jString::string_replace(string strold, const char c, const char r)
{
    int N = strold.length();
    char* str = (char *)strold.data();
    char* tmp = new char[N];
    int i=0;
    while(str[i])
    {
        tmp[i] = str[i]; if(str[i]==c) tmp[i]=r;
        i++;
    }
    string u(tmp,N);
    delete [] tmp;
    return u;
}

string jString::string_replace(string strold, const string subold, const string subnew)
{
    string strnew = strold;
    int loc1 = strold.find(subold);
    int loc2 = loc1 + subold.length();
    while(loc1 >= 0 && loc2 < string::npos)
    {
        string tmp1 = strnew.substr(0,loc1);
        string tmp2 = strnew.substr(loc2,strnew.length() - loc2);
        strnew = tmp1+subnew+tmp2;
        loc1 = strnew.find(subold);
        loc2 = loc1 + subold.length();
    }
    return strnew;
}

string jString::string_trim(string s)
{
    int N = s.length();
    char* str = (char *)s.data();
    char* tmp = new char[N];
    int i=0;
    int j=0;
    for(i=0;i<N;i++)
    {
        if(str[i] != ' ') tmp[j++] = str[i];
    }
    string u(tmp,j);
    delete [] tmp;
    return u;
}

string jString::string_whiteOut(string s, const string r)
{
    int N = s.length();
    int M = r.length();
    char* str = (char *)s.data();
    char* rem = (char *)r.data();
    char* tmp = new char[N];
    int i=0;
    while(str[i])
    {
        tmp[i] = str[i];
        for(int j=0;j<M;j++) if(str[i]==r[j]) tmp[i]=' ';
        i++;
    }
    string u(tmp,N);
    delete [] tmp;
    return u;
}

string jString::string_clean(const string strin)
{
    return string_clean(strin," ,;:=#$\t\n\r\f");
}

string jString::string_clean(const string strin, const string separ)
{
    string cleaned = strin;
    int beg, loc1, loc2 = 0;
    bool comment = false;
    if(cleaned.compare(0,1,"#")==0) comment = true;
    if(cleaned.length() > 1)
    {
        int beg = cleaned.find_first_not_of(' ',0);
        if(cleaned.compare(beg,2,"//")==0) comment = true;
    }
    if(comment) { return ""; }
    // remove trailing comments
    string ss = "//";
    beg = cleaned.find(ss);
    //beg = cleaned.find_first_of('//',0);
    if(beg > 0 && beg != string::npos)
    {
        if(cleaned.compare(beg,2,"//")==0) cleaned.erase(beg,cleaned.length());
    }
    beg = cleaned.find_first_of('#',0);
    if(beg > 0) cleaned.erase(beg,cleaned.length());
    // remove leading and trailing separators
    loc1 = cleaned.find_first_not_of(separ);
    loc2 = cleaned.find_last_not_of( separ);
    cleaned = cleaned.substr(loc1,loc2-loc1+1);
    return cleaned;
}

bool jString::contains(const string s, const string match)
{
    bool ret = false;
    if(s.find(match) < string::npos) ret=  true;
    return ret;
}

bool jString::endsWith(const string s, const string match)
{
    int beg =  s.length()-match.length(); 
    string sub = s.substr(beg,s.length());
    return (s.substr(beg,s.length())==match);    
}

string jString::format(const string f, int value)
{
    int beg = f.find_first_of('%',0);
    int end = f.find_first_of('d',beg);
    string s = f;
    if(beg >= 0 && end < string::npos)
    {
        string tmp1 = f.substr(beg,end-beg+1);
        string tmpw = f.substr(beg+1,end-beg);
        int width = atoi(tmpw.c_str());
        string tmp2 = jString::toString(value,width);
        s = jString::string_replace(s, tmp1, tmp2);
    }
    return s;
}

string jString::format(const string f, float value)
{
    int beg = f.find_first_of('%',0);
    int dot = f.find_first_of('.',beg);
    int end = f.find_first_of('f',beg);
    string s = f;
    if(beg >= 0 && end < string::npos)
    {
        if(dot<=beg) dot=end;
        string tmp1 = f.substr(beg,end-beg+1);
        string tmpw = f.substr(beg+1,dot-beg);
        string tmpp = f.substr(dot+1,end-dot);
        int prec = atoi(tmpp.c_str());
        string tmp2 = jString::toString(value,prec);
        s = jString::string_replace(s, tmp1, tmp2);
    }
    return s;
}

string jString::format(const string f, double value)
{
    int beg = f.find_first_of('%',0);
    int dot = f.find_first_of('.',beg);
    int end = f.find_first_of('e',beg);
    string s = f;
    if(beg >= 0 && end < string::npos)
    {
        if(dot<=beg) dot=end;
        string tmp1 = f.substr(beg,end-beg+1);
        string tmpw = f.substr(beg+1,dot-beg);
        string tmpp = f.substr(dot+1,end-dot);
        int prec = atoi(tmpp.c_str());
        string tmp2 = jString::toString(value,prec);
        s = jString::string_replace(s, tmp1, tmp2);
    }
    return s;
}

string jString::format(const string f, string value)
{
    int beg = f.find_first_of('%',0);
    int end = f.find_first_of('s',beg);
    string s = f;
    if(beg >= 0 && end < string::npos)
    {
        string tmp1 = f.substr(beg,end-beg+1);
        string tmpw = f.substr(beg+1,end-beg);
        int width = atoi(tmpw.c_str());
        ostringstream outstr;
        outstr.setf(ios_base::right, ios_base::adjustfield);
        outstr.width(width);
        outstr << value;
        string tmp2 = outstr.str();
        s = jString::string_replace(s, tmp1, tmp2);
    }
    return s;
}

    
Parsed::Parsed() {field = NULL; ftype = NULL; value = NULL; };
Parsed::Parsed(string sin)
{
	field = NULL; ftype = NULL; value = NULL;
    string sep =" ,;:=#$\t\n\r\f";
    init(sin,sep,true); 
}
Parsed::Parsed(string sin, bool par)
{
	field = NULL; ftype = NULL; value = NULL;
    string sep =" ,;:=#$\t\n\r\f";
    init(sin,sep,par); 
}
Parsed::Parsed(string sin, string sep)
{
	field = NULL; ftype = NULL; value = NULL;
    init(sin,sep,true);  
}
Parsed::Parsed(string sin, string sep, bool par)
{
	field = NULL; ftype = NULL; value = NULL;
    init(sin,sep,par);  
}
Parsed::Parsed(jString js)
{
	field = NULL; ftype = NULL; value = NULL;
    string sep =" ,;:=#$\t\n\r\f";
    init(js.str(),sep,true); 
}
Parsed::Parsed(jString js, bool par)
{
	field = NULL; ftype = NULL; value = NULL;
    string sep =" ,;:=#$\t\n\r\f";
    init(js.str(),sep,par); 
}
Parsed::Parsed(jString js, string sep)
{
	field = NULL; ftype = NULL; value = NULL;
    init(js.str(),sep,true);  
}
Parsed::Parsed(jString js, string sep, bool par)
{
	field = NULL; ftype = NULL; value = NULL;
    init(js.str(),sep,par);  
}

void Parsed::init(string sin, string sep, bool par)
{
    ignore_parenthesis = par;
    
    if(ignore_parenthesis)
    {
        strin = sin;
    }
    else
    {
        char* str = (char *)sin.data();
        vector<char> tmp;
        int ics = 0;
        int i=0;
        while(str[i])
        {
            switch (str[i])
            {                
                case '(':
                    tmp.push_back(' '); tmp.push_back(str[i]); ics=1; break;
                case ')':
                    tmp.push_back(' '); tmp.push_back(str[i]); ics=0; break;
                case ',':
                    if(ics>0) { tmp.push_back(')'); tmp.push_back(' '); tmp.push_back('('); }
                    else tmp.push_back(str[i]);
                    break;
                case ' ':
                    if(ics == 0) tmp.push_back(' ');
                    break;
                default:
                    tmp.push_back(str[i]);                    
            }
            i++;
        }
        string strin(tmp.data(),tmp.size());
    }
    
    nflds = 0;
    type = "COMM";
    separ = sep;
    length= strin.length();
    
    vector<string> tokens = jString::string_parse(strin,separ);
    nflds = tokens.size();
    field = new jString[nflds];
    ftype = new jString[nflds];
    value = new  double[nflds];
    for(int n=0;n<nflds;n++)
    {
        field[n] = jString(tokens[n]);
        ftype[n] = jString("NUMB");
        try
        {
            value[n] = boost::lexical_cast<double>(field[n].chr());
        }
        catch(std::exception &)
        {
            ftype[n] = jString("WORD");
            value[n] = 0.0;
        }
    }
    
    
}

vector<jString> Parsed::read_var(string name)
{
	vector<jString> vals;
	if (!field[0].equals(name)) { printf("found %s instead of %s\n",field[0].str().c_str(),name.c_str()); abort(); }
	if (nflds==1) { printf("do not find any value of %s\n",name.c_str()); abort();}
	for(int i=1;i<nflds;i++)
		vals.push_back(field[i]);
	return vals;
}

void Scan::open()
{
    fi.open(filename.c_str(),ios_base::in);
    if(fi.fail()) { cout << "file:" << filename <<": could not be opened" << endl;  abort(); }
    fi.seekg(0,ios_base::beg);
}
void Scan::open(string path)
{
    filename = path; 
    open();
}
void Scan::close()
{
    fi.close();
}

jString Scan::readLine()
{
    jString jline;
    if(!fi.eof())
    {
        string line;
        if(getline(fi,line) ) jline = jString(line);
    }
    return jline;
}

jString Scan::nextLine()
{
    jString line;
    if(next_line < input.size())  line = input[next_line++]; 
    return line;
}

jString Scan::Line(int i)
{
    return input[i];
}

jString Scan::Line(string key)
{
	rewind();
	while(!end_of_scan())
	{
		jString jline = nextLine().clean();
		Parsed parsed(jline);
		if (parsed.field[0].equals(key))
			return jline;
	}
	printf("cannot find line starting with %s\n",key.c_str());
    return jString("");
}

void Scan::iterNext() { next_line++; }
void Scan::iterBack() { next_line--; }

void Scan::readInputList()
{
    // check file exists and position at start
    fi.open(filename.c_str(),ios_base::in);
    if(fi.fail()) { cout << "file:" << filename <<": could not be opened" << endl;  abort(); }
    fi.seekg(0,ios_base::beg);

    string line;
    // read all lines - ignore comments
    while(!fi.eof())
    {
        if(getline(fi,line) )
        {            
            jString jline = jString(line);
            jline = jline.clean().whiteOut("()").trim();
            if(jline.length()==0)  continue;
            
            input.push_back(jline);

        } // read next line
    }
    fi.close();    
    next_line = 0;
}

void Scan::readInputList(string def_KW)
{
    // check file exists and position at start
    fi.open(filename.c_str(),ios_base::in);
    if(fi.fail()) { cout << "file:" << filename <<": could not be opened" << endl;  abort(); }
    fi.seekg(0,ios_base::beg);

    bool def = false;
    string line;

    // read all lines and select those for correct definition - ignore comments
    while(!fi.eof())
    {
        if(getline(fi,line) )
        {            
            jString jline = jString(line);
            jline = jline.clean().whiteOut("()").trim();
            if(jline.length()==0)  continue;
            
            vector<jString> tokens = jline.toUpperCase().parse();
            if(tokens.size() == 0) continue;
            
            int nf = 0;
            while(nf < tokens.size())
            {
                jString tok = tokens[nf];
                // end definitions
                if( tok.contains("END") )
                {  
                    if(def) 
                    {
                        jString typ = tokens[++nf];
                        if(typ.contains("DEF") || typ.contains(def_KW)) def = false;
                    }
                    break; 
                }
                
                // new definition
                if( tok.contains("DEF") )
                {
                    jString typ = tokens[++nf];     
                    if(typ.contains(def_KW)) def = true;
                }
                nf++;
            }
            tokens.clear();
            if(def) input.push_back(jline);

        } // read next line
    }
    fi.close();    
    next_line = 0;
}

// formatting utilities
using boost::format;
string myiolib::myformat(string a)
{
    int n = a.length();
    string s = boost::io::str(format("%d") % n);
    string form = "%" + s + "s";
    s = boost::io::str(format(form) % a);
    return s;
}
string myiolib::myformat(int a)
{
    int n = 1; if(abs(a) > 0) n = (int)log10(1.*abs(a));
    if(a < 0) n++;
    string s = boost::io::str(format("%d") % n);
    string form = "%" + s + "d";
    s = boost::io::str(format(form) % a);
    return s;
}
string myiolib::myformat(string a, int n)
{
    string s = boost::io::str(format("%d") % n);
    string form = "%" + s + "d";
    s = boost::io::str(format(form) % a);
    return s;
}


