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

int main()
{
  const string dataFolder("doc_data");
  ifstream ifs(dataFolder+"/DocList.txt");
  //intro(), website(");
  vector<string> fcnList, introList, websiteList;
  string line;
  while (ifs >> line) {
    fcnList.push_back(line);
  }
  ifs.close();
  ifs.open(dataFolder+"/intro.m");
  while (getline(ifs, line)) {
    introList.push_back(line);
  }
  ifs.close();
  ifs.open(dataFolder+"/website.m");
  while (getline(ifs, line)) {
    websiteList.push_back(line);
  }
  ifs.close();
  ofstream gail("GAIL_t.m");
  for (const auto &s : introList) {
    gail << s << "\n";
  }
  gail << "%% Functions\n";
  gail << "%\n";
  gail << "% <html>" << endl;
  for (const auto &s : fcnList) {
    gail << "% <a href=\"help_" << s << ".html\">" << s << "</a>\n";
  }
  gail << "% </html>\n";
  gail << "%\n";
  gail << "%" << endl;
  for (const auto &s : websiteList) {
    gail << s << "\n";
  }
  gail << std::flush;
  gail.close();
}
