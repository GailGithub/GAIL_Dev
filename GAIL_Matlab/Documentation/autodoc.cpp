/*
 * This program automatically generates .m files that can be published to html documentation.
 * Refer to https://sites.google.com/site/gailteam1/project-updates/howtouseautodoctocreatehtmldocumentation
 * for details.
*/ 

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::endl;
using std::flush;
using std::find;

string upperString(const string &) noexcept;
string lowerString(const string &) noexcept;

int main()
{
  const string dataFolder("doc_data");
  ifstream ifs(dataFolder + "/DocList.txt");
  vector<string> fcnList, fcnName, uFcnName, introList, websiteList, fcnDoc;
  string line, word;
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
  ofstream gail("GAIL_raw.m"), helptoc("html/helptoc_raw.xml"), funclist("funclist_raw.m"), ofs;
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
    if (!s.empty() && *(s.cend() - 1) != ':') {
      fcnName.push_back(s);
      uFcnName.push_back(upperString(s));
    }
  }
  for (const auto &s : fcnList) {
    if (s.empty()) {
      funclist << "%\n";
    } else if (*(s.cend() - 1) == ':') {
      funclist << "%% " << s.substr(0,s.size() - 1) << "\n%\n";
    } else {
      helptoc << "            <tocitem target=\"help_" << s << ".html\">" << s << "</tocitem>\n";
      gail << "% <a href=\"help_" << s << ".html\">" << s << "</a>\n";
      funclist << "% <html>\n% <a href=\"help_" << s << ".html\">" << s << "</a>\n% </html>\n%\n";
      ifs.open("../Algorithms/" + s + ".m");
      while (getline(ifs, line) && line != "") {
	fcnDoc.push_back(line);
      }
      ifs.close();
      ofs.open("help_" + s + "_raw.m");
      ofs << "%% " << s << "\n% ";
      auto space1 = find(fcnDoc[1].cbegin(), fcnDoc[1].cend(), ' ');
      fcnDoc[1] = fcnDoc[1].substr(space1 - fcnDoc[1].cbegin() + 1, fcnDoc[1].size());
      ofs << fcnDoc[1];
      auto emptyLine1 = find(fcnDoc.cbegin(), fcnDoc.cend(), "%");
      for (auto iter = fcnDoc.begin() + 2; iter != emptyLine1; ++iter) {
	if (iter != emptyLine1 -1) {
	  ofs << "\n% " << (*iter).substr(1);
	} else {
	  ofs << "\n% " << (*iter).substr(1);
	}
      }
      ofs << "\n%% Syntax" << endl;
      auto inputArg = find(++emptyLine1, fcnDoc.cend(), "%   Input Arguments");
      us = upperString(s);
      {
	decltype(fcnDoc.size()) cnt = 0;
	for (auto iter = emptyLine1; iter != inputArg; ++iter) {
	  auto lPos = (*iter).find(" = " + us + "(");
	  if (lPos != string::npos) {
	    ++cnt;
	    string sReplace = *iter;
	    sReplace.replace(lPos + 3, s.size(), "*" + s + "*");
	    auto rPos = sReplace.find_first_of(')', lPos);
	    if (cnt == 1) {
	      ofs << "% " << sReplace.substr(4, rPos - 3) << "\n";
	    } else {
	      ofs << "%\n% " << sReplace.substr(4, rPos - 3) << "\n";
	    }
	  }
	}
      }
      ofs << "%% Description\n%";
      for (auto iter = emptyLine1; iter != inputArg; ++iter) {
	auto lPos = (*iter).find(" = " + us + "(");
	if (lPos == string::npos) {
	  if ((*iter).size() > 4) {
	    ofs << "\n%  " <<  (*iter).substr(4);
	  } else {
	    ofs << "\n" << *iter;
	  }
	} else {
	  string sReplace = *iter;
	  sReplace.replace(lPos + 3, s.size(), "*" + s + "*");
	  ofs << "\n% " << sReplace.substr(4);
	}
      }

      ofs << "\n% *Input Arguments*\n%" << endl;
      auto optInputArg = find(++inputArg, fcnDoc.cend(), "%   Optional Input Arguments");
      auto outputArg = find(inputArg, fcnDoc.cend(), "%   Output Arguments");
      if (optInputArg != fcnDoc.cend()) {
	for (auto iter = ++inputArg; iter != optInputArg; ++iter) {
	  if (iter->size() > 6) {
	    if (iter->find(" --- ") != string::npos) {
	      ofs << "% * " << iter->substr(6);
	    } else {
	      ofs << "\n%  " << iter->substr(6);
	    }
	  } else {
	    ofs << "\n" << *iter << "\n";
	  }
	}
	ofs << "% *Optional Input Arguments*\n%" << endl;
	for (auto iter = optInputArg + 2; iter != outputArg; ++iter) {
	  if (iter->size() > 6) {
	    if (iter->find(" --- ") != string::npos) {
	      ofs << "% * " << iter->substr(6);
	    } else {
	      ofs << "\n%  " << iter->substr(6);
	    }
	  } else {
	    ofs << "\n" << *iter << "\n";
	  }
	}
      } else {
	for (auto iter = ++inputArg; iter != outputArg; ++iter) {
	  if (iter->size() > 6) {
	    if (iter->find(" --- ") != string::npos) {
	      ofs << "% * " << iter->substr(6);
	    } else {
	      ofs << "\n%  " << iter->substr(6);
	    }
	  } else {
	    ofs << "\n" << *iter << "\n";
	  }
	}
      }

      ofs << "% *Output Arguments*\n%" << endl;
      auto guarantee = find(++outputArg, fcnDoc.cend(), "%  Guarantee");
      for (auto iter = ++outputArg; iter != guarantee; ++iter) {
	if ((*iter).size() > 6) {
	  auto pos = (*iter).find(" --- ");
	  if (pos != string::npos) {
	    ofs << "% * " << iter->substr(6);
	  } else {
	    ofs << "\n%  " << iter->substr(6);
	  }
	} else {
	  ofs << "\n" << *iter << "\n";
	}
      }
      ifstream fcnData(dataFolder + "/" + s + "_data.m");
      while (getline(fcnData, line)) {
	ofs << line << "\n";
      }
      ofs << "%% See Also\n%" << endl;
        


      auto see = find_if(++guarantee, fcnDoc.cend(), [](const string &a) { return a.size() >= 12 && a.substr(4,8) == "See also"; });
      istringstream sa((*see).substr(13));
      while (sa >> word) {
        if (*(word.end() - 1) == ',') {
	  word.erase(word.end() - 1);
	}
	auto num = find(uFcnName.cbegin(), uFcnName.cend(), word);
    string url_pre = "\"help_";
    if (std::find(fcnList.begin(), fcnList.end(), fcnName[num - uFcnName.begin()]) != fcnList.end()) {
      url_pre = "\"http://www.mathworks.com/help/matlab/ref/";
    }
    ofs << "% <html>\n% <a href=" + url_pre << fcnName[num - uFcnName.begin()] << ".html\">" << fcnName[num - uFcnName.begin()] << "</a>\n% </html>\n%\n";

      }
      ofs << "%% References" << endl;
      auto ref = find(++see, fcnDoc.cend(), "%  References");
      for (auto iter = ++ref; iter != fcnDoc.cend(); ++iter) {
	if ((*iter).size() > 4) {
	  ofs << "% " << (*iter).substr(4) << "\n";
	} else {
	ofs << *iter << "\n";
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
  std::cout << "autodoc: Automatic documentation is completed." << endl;
}

string upperString(const string &s) noexcept
{
  string uStr(s.size(),' ');
  for (string::size_type i = 0;i != s.size(); ++i) {
    uStr[i] = toupper(s[i]);
  }
  return uStr;
}

string lowerString(const string &s) noexcept
{
  string lStr(s.size(),' ');
  for (string::size_type i = 0;i != s.size(); ++i) {
    lStr[i] = tolower(s[i]);
  }
  return lStr;
}
