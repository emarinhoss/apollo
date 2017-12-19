/*
 * File:   myIO.h
 * Author: cambier
 *
 * Created on December 6, 2011, 2:08 PM
 */

#ifndef MYIO_H
#define	MYIO_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <vector>
#include <map>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#define BUF_SIZE  1024

using namespace std;

namespace myiolib
{

    /*
     * IO utility
     */
    class JIO
    {        
    public:
        FILE* f;
        bool opened;

        JIO(): f(NULL),opened(false) { }
       ~JIO() { if(f != NULL)  {  if(opened) fclose(f); f = NULL; } opened = false; }

        void open(string filename, string type);
        void open(const char* filename, const char* type);
        void close();
        bool exists(string);

        void writeInt   (int*    pi);
        void writeFloat (float*  pr);
        void writeDouble(double* pr);
        void writeBool  (bool*   pb);
        void writeString(string s);
        void writeInt   (int    i);
        void writeFloat (float  r);
        void writeDouble(double r);
        void writeBool  (bool   b);

        int    readInt()   ;
        float  readFloat() ;
        double readDouble();
        bool   readBool()  ;
        string readString();
        void writeIntArray   (int N, int*    pa);
        void writeFloatArray (int N, float*  pa);
        void writeDoubleArray(int N, double* pa);
        int*   readIntArray    (int N);
        float* readFloatArray  (int N);
        double* readDoubleArray(int N);

        void dump();
        static int makeDirectory(string sdir);
        static int checkDirectory(string sdir);
        static bool is_Directory(string sdir);
        static bool is_File(string sdir);

    };



    class CppIOinp
    {
    private:
        ifstream fi;
    public:
        CppIOinp() {;}
        CppIOinp(string filename) {  fi.open(filename.c_str(),ios_base::in  | ios_base::binary);  }
       ~CppIOinp() {  fi.close(); }
        float  readFloat () {  float  f1; fi.read((char *) &f1,sizeof(float )); return f1;  }
        double readDouble() {  double d1; fi.read((char *) &d1,sizeof(double)); return d1;  }
        int    readInt   () {  int i1;    fi.read((char *) &i1,sizeof(int))   ; return i1;  }
        string readString()
        {
            char name[512];
            int i=0;
            char c; fi.get(c); while (c) {  name[i++] = c; fi.get(c); }
            name[i]=c; string s(name);
            return s;
        }

    };

    class CppIOout
    {
    private:
        ofstream fo;
    public:
        CppIOout() {;}
        CppIOout(string filename) {  fo.open(filename.c_str(),ios_base::out | ios_base::binary);  }
       ~CppIOout() {  fo.close(); }
        inline void writeDouble(double d1)  {  fo.write((char *) &d1,sizeof(double));  }
        inline void writeFloat (float  f1)  {  fo.write((char *) &f1,sizeof(float ));  }
        inline void writeInt   (int    i1)  {  fo.write((char *) &i1,sizeof(int   ));  }
        inline void writeString(string s)   {  fo.write(s.c_str(),s.length()); fo.put('\0');  }
    };

    /**
     * Java-like string object
     */
    class jString
    {
    private:
        string s;
    public:
        jString(): s("")  { }
        jString(string s) {this->s = s; }
       ~jString() { }

        void reset() { s=""; }
        void reset(char* a) { int i=0; while(a[i]) { this->s[i] = a[i]; i++; } }
        string str() { return s; }
        char*  chr() { return (char*) s.c_str(); }
        char*  chr(int is, int ie) { return (char*) s.substr(is,ie).c_str(); }
        int length() { return (int)s.length(); }
        inline void operator=(const string sr) { s=sr; }
        inline bool operator==(jString c) { return (s==c.str()); }
        inline bool operator==( string c) { return (s==c); }

        jString trim();
        jString clean();
        jString clean(const string sep);
        vector<jString> parse();
        vector<jString> parse(const string sep);
        jString toUpperCase();
        jString whiteOut(const string);
        jString concat(jString);
        jString concat (string);
        jString substring(int);
        jString substring(int,int);
        bool contains(string);
        bool matches(int, string);
        bool equals (jString);
        bool equals ( string);
        void replace(const char  ,const char);
        void replace(const string,const string);
        int indexOf(const char);
        int indexOf(const string);
        int lastIndexOf (const char);
        int lastIndexOf (const string);
        bool startsWith(string);
        bool endsWith(string);
        char charAt(int);
        void print() { printf("%s\n",s.c_str()); };
        
        int    iValue();
        double dValue();
        float  fValue();
        double rValue();

        static string toString(int value);
        static string toString(int value, int width);
        static string toString(int value, int width, char filler);
        static string toString(float  value);
        static string toString(float  value, int prec);
        static string toString(double value);
        static string toString(double value, int prec);
        static vector<string> string_parse(const string sin);
        static vector<string> string_parse(const string sin, const string sep);
        static string string_replace(string strold, const string subold, const string subnew);
        static string string_replace(string strold, const char c, const char r);
        static string string_trim(string strin);
        static string string_clean(const string strin);
        static string string_clean(const string strin, const string separ);
        static string string_whiteOut(string s, const string r);
        static bool contains(const string s, const string match);
        static bool endsWith(const string s, const string match);
        
        static string format(const string,int);
        static string format(const string,float);
        static string format(const string,double);
        static string format(const string,string);

    }; //end String class


    
    class Parsed
    {
    public:
        int length,nflds;
        string type,strin,separ;
        jString* field;
        jString* ftype;
        double*  value;
        bool ignore_parenthesis;
        
        ~Parsed() { if (field) delete [] field; if (ftype) delete [] ftype; if (value) delete [] value; }
        Parsed();
        Parsed(string);
        Parsed(string,string);
        Parsed(string,bool);
        Parsed(string,string,bool);
        Parsed(jString);
        Parsed(jString,string);
        Parsed(jString,bool);
        Parsed(jString,string,bool);
        
        void init(string,string,bool);

        vector<jString> read_var(string);
    };
    

    /*
     * System-information (equivalent of Java jFileSys)
     */
    class jFileSys
    {
    private:
    	// Previously incapable of compiling due to private 'init()'
//        void init();
//        void kill();
//        void init(string);

    public:
        void init();
        void kill();
        void init(string);

        JIO jio;
        string data_dir;  // data directory (read-only)
        string curr_dir;  // current directory
        string home_dir;  // home directory
        string brun_dir;  // binary output data
        string case_dir;  // case directory( where input and output data is)
        string username;  // user name
        string casename;  // case-name
        vector<string> data_subd; // sub-directories in Database
        
        string thrmfile;  // binary thermo-chemical file

        jFileSys();
       ~jFileSys(); 
        jFileSys(string name);
        
        void display();

    };   
    

    class Scan
    {
    private:
        vector<string>  def_name;
        vector<string>* val_name;
        int NvalMax;
        bool nodeData;
        bool cellData;
        bool faceData;

    public:
        jFileSys jsys;
        ifstream fi;
        string filename;
        vector<jString> input;
        int next_line;

        string where;
        int NdataSetup;
        vector<int>     pointsSet;
        vector<int*>    selectSet;
        vector<float*>  fvalueSet;
        vector<string*> svalueSet;

        Scan() { }
        Scan(const Scan& orig) {;}
       ~Scan() { close(); }
        Scan(string path)
        {
            filename = path;
            iocheck();
        }

        bool end_of_scan() { if (next_line == (int) input.size()) return true; else return false; };
        void iocheck()
        {
            fi.open(filename.c_str(),ios_base::in);
            if(fi.fail()) { cout << "input file:" << filename <<": could not be opened" << endl; abort(); }
            fi.close();
        }
        void readInputList();
        void readInputList(string);
        jString readLine();
        jString nextLine();
        jString Line(int);
        jString Line(string);
        void iterNext();
        void iterBack();
        void rewind() { while (next_line != 0) iterBack(); }
        void open(string);
        void open();
        void close();

    };


    //formatting utilities

    string myformat(string a);
    string myformat(string a, int n);
    string myformat(int a);
    string myformat(int a, int n);

}
    
#endif	/* MYIO_H */



