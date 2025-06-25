#include "printStars.h"
#include <Gaudi/Algorithm.h>
#include <fmt/core.h>
#include <string>

void printInStars(const Gaudi::Algorithm* thisAlg, const std::string& msg, const int lineWidth) {
  thisAlg->debug() << fmt::format("{:*^{}}", "", lineWidth) << endmsg;
  thisAlg->debug() << fmt::format("{:*^{}}", msg, lineWidth) << endmsg;
  thisAlg->debug() << fmt::format("{:*^{}}", "", lineWidth) << endmsg;
}
