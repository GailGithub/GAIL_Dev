#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
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
  while (ifs >> line) {
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
  ofstream gail("GAIL_t.m"), ofs;
  for (const auto &s : introList) {
    gail << s << "\n";
  }
  gail << "%% Functions\n" << "%\n" << "% <html>" << endl;
  for (const auto &s : fcnList) {
    gail << "% <a href=\"help_" << s << ".html\">" << s << "</a>\n";
    ifs.open("../Algorithms/" + s + ".m");
    while (getline(ifs, line) && line != "") {
      // string::size_type pos = 0;
      // auto uStr = upperString(s);
      // while ((pos = line.find(uStr, pos)) != string::npos) {
      // 	line.replace(pos, s.size(), uStr);
      // 	pos += uStr.length();
      // }
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
    ofs << "|\n%% Syntax" << endl;
    ofs.close();
    fcnDoc.clear();
  }
  gail << "% </html>\n" << "%\n" << "%" << endl;
  for (const auto &s : websiteList) {
    gail << s << "\n";
  }
  gail << flush;
  gail.close();
}

string upperString(const string &s) noexcept
{
  string uStr;
  for (string::size_type i = 0;i != s.size();++i) {
    uStr[i] = toupper(s[i]);
  }
  return uStr;
}
