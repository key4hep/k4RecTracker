#ifndef TESTFWCORE_HELLOWORLDALG
#define TESTFWCORE_HELLOWORLDALG

#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

class HelloWorldAlg : public GaudiAlgorithm {
public:
  explicit HelloWorldAlg(const std::string&, ISvcLocator*);
  virtual ~HelloWorldAlg();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // member variable
  Gaudi::Property<std::string> theMessage{this, "PerEventPrintMessage", "Hello ",
                                          "The message to printed for each Event"};
};

#endif /* TESTFWCORE_HELLOWORLDALG */
