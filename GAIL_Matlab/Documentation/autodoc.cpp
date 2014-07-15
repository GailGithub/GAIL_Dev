#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <cctype>
//#include <typeinfo>
#include <algorithm> 
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
//using std::fstream;
using std::cout;
using std::endl;
using std::flush;
using std::find;

string upperString(const string &) noexcept;

int main()
{
  const string dataFolder("doc_data");
  ifstream ifs(dataFolder + "/DocList.txt");
  vector<string> fcnList, introList, websiteList, fcnDoc;
  string line;
  while (getline(ifs,line)) {
    fcnList.push_back(line);
  }
  ifs.close();
  ifs.open(dataFolder + "/intro.m");
  while (getline(ifs, line)) {
    introList.push_back(line);
  }
  ifs.close();
  ifs.open(dataFolder + "/website.m");
  while (getline(ifs, line)) {
    websiteList.push_back(line);
  }
  ifs.close();
  ofstream gail("GAIL_t.m"), helptoc("html/helptoc_t.xml"), funclist("funclist_t.m"), ofs;
  helptoc << "<?xml version='1.0' encoding='ISO-8859-1' ?>\n\n<toc version=\"1.0\">\n\n"
          << "<tocitem target=\"GAIL.html\">GAIL Toolbox\n"
          << "    <tocitem target=\"funclist.html\" image=\"HelpIcon.FUNCTION\">Functions" << endl;
  for (const auto &s : introList) {
    gail << s << "\n";
  }
  gail << "%% Functions\n" << "%\n" << "% <html>" << endl;
  funclist << "%% Functions\n%" << endl;
  string us;
  for (const auto &s : fcnList) {
    if (s.empty()) {
      funclist << "%\n";
    } else if (*(s.cend()-1) == ':') {
      funclist << "%% " << s.substr(0,s.size() - 1) << "\n%\n";
    } else {
      helptoc << "            <tocitem target=\"help_" << s << ".html\">" << s << "</tocitem>\n";
      gail << "% <a href=\"help_" << s << ".html\">" << s << "</a>\n";
      funclist << "% <html>\n% <a href=\"help_" << s << ".html\">meanMC_g</a>\n% </html>\n%\n";
      ifs.open("../Algorithms/" + s + ".m");
      while (getline(ifs, line) && line != "") {
	fcnDoc.push_back(line);
      }
      ifs.close();
      ofs.open("help_" + s + ".m");
      ofs << "%% " << s << "\n% |";
      auto space1 = find(fcnDoc[1].cbegin(), fcnDoc[1].cend(), ' ');
      fcnDoc[1] = fcnDoc[1].substr(space1 - fcnDoc[1].cbegin() + 1, fcnDoc[1].size());
      auto emptyLine1 = find(fcnDoc.cbegin(), fcnDoc.cend(), "%");
      for (auto iter = ++fcnDoc.begin(); iter != emptyLine1; ++iter) {
	if (iter != emptyLine1 -1) { 
	  ofs << *iter << "\n";
	} else {
	  ofs << *iter;
	}
      }
      ofs << ".|\n%% Syntax" << endl;
      auto guarantee = find(++emptyLine1, fcnDoc.cend(), "%  Guarantee");
      us = upperString(s);
      {
	decltype(fcnDoc.size()) cnt = 0;
	for (auto iter = emptyLine1; iter != guarantee; ++iter) {
	  auto lPos = (*iter).find(" = " + us + "(");
	  if (lPos != string::npos) {
	    ++cnt;
	    auto rPos = (*iter).find_first_of(')', lPos);
	    if (cnt == 1) {
	      ofs << "% " << (*iter).substr(4, rPos - 3) << "\n";
	    } else {
	      ofs << "%\n% " << (*iter).substr(4, rPos - 3) << "\n";
	    }
	  }
	}
      }
      ofs.flush();
      ofs.close();
      fcnDoc.clear();
    }
  }
  helptoc << "        </tocitem>\n    </tocitem>\n</toc>" << endl;
  gail << "% </html>\n" << "%\n" << "%" << endl;
  for (const auto &s : websiteList) {
    gail << s << "\n";
  }
  gail << flush;
  funclist << flush;
  helptoc.close();
  gail.close();
  funclist.close();
}

string upperString(const string &s) noexcept
{
  string uStr(s.size(),' ');
  for (string::size_type i = 0;i != s.size(); ++i) {
    uStr[i] = toupper(s[i]);
  }
  return uStr;
}
