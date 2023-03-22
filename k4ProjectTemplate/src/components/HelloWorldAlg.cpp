#include "HelloWorldAlg.h"
#include "GaudiKernel/MsgStream.h"

DECLARE_COMPONENT(HelloWorldAlg)

HelloWorldAlg::HelloWorldAlg(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

HelloWorldAlg::~HelloWorldAlg() {}

StatusCode HelloWorldAlg::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode HelloWorldAlg::execute() {
  info() << endmsg;
  info() << endmsg;
  info() << theMessage << endmsg;
  info() << endmsg;
  info() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode HelloWorldAlg::finalize() { return GaudiAlgorithm::finalize(); }
