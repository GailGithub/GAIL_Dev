#include <iostream>
#include <fstream>
#include <string>
#include <list>
using std::list;
using std::string;
using std::ifstream;
using std::ofstream;
//using std::fstream;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  const string dataFolder("doc_data");
  ifstream doclist(dataFolder+"/DocList.txt"), intro(dataFolder+"/intro.m"), website(dataFolder+"/website.m");
  list<string> fcnList, introList, websiteList;
  string line;
  while (doclist >> line) {
    fcnList.push_back(line);
  }
  doclist.close();
  while (getline(intro, line)) {
    introList.push_back(line);
  }
  intro.close();
  while (getline(website, line)) {
    websiteList.push_back(line);
  }
  website.close();
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
  gail << "%\n";
  gail << "%" << endl;
  for (const auto &s : websiteList) {
    gail << s << "\n";
  }
  gail << std::flush;
}
