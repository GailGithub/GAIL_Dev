#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
//using std::fstream;
using std::cout;
using std::endl;
using std::flush;

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
  gail << "%% Functions\n";
  gail << "%\n";
  gail << "% <html>" << endl;
  for (const auto &s : fcnList) {
    gail << "% <a href=\"help_" << s << ".html\">" << s << "</a>\n";
    ifs.open("../Algorithms/" + s + ".m");
    while (getline(ifs, line) && line != "") {
      fcnDoc.push_back(line);
    }
    cout << fcnDoc.size() << endl;
    fcnDoc.clear();
    ifs.close();
    ofs.open("help_" + s + ".m");
    ofs << "%% " << s << "\n";
    ofs << flush;
    ofs.close();
  }
  gail << "% </html>\n";
  gail << "%\n";
  gail << "%" << endl;
  for (const auto &s : websiteList) {
    gail << s << "\n";
  }
  gail << flush;
  gail.close();
}
