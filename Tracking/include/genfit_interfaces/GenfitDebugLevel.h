#ifndef GENFIT_DEBUG_LEVEL_H
#define GENFIT_DEBUG_LEVEL_H

#include <cstring>
#include <streambuf>

#include <IO.h>
#include <TError.h>
#include <TGeoManager.h>

namespace {

class NullStreamBuffer final : public std::streambuf {
protected:
  int_type overflow(int_type character) override { return traits_type::not_eof(character); }
};

ErrorHandlerFunc_t previousRootErrorHandler = nullptr;

void filteredRootErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
  const bool isNonPositiveDefiniteCholesky = location != nullptr && message != nullptr &&
                                             std::strstr(location, "TDecompChol::Decompose") != nullptr &&
                                             std::strstr(message, "matrix not positive definite") != nullptr;

  if (isNonPositiveDefiniteCholesky) {
    return;
  }

  if (previousRootErrorHandler != nullptr) {
    previousRootErrorHandler(level, abort, location, message);
  } else {
    DefaultErrorHandler(level, abort, location, message);
  }
}

class ScopedFitDiagnostics final {
public:
  explicit ScopedFitDiagnostics(bool suppress) : m_suppress(suppress) {
    if (!m_suppress) {
      return;
    }

    m_genfitErrorBuffer = genfit::errorOut.rdbuf(&m_sink);
    m_genfitDebugBuffer = genfit::debugOut.rdbuf(&m_sink);
    m_previousRootErrorHandler = SetErrorHandler(filteredRootErrorHandler);
    previousRootErrorHandler = m_previousRootErrorHandler;
  }

  ~ScopedFitDiagnostics() {
    if (!m_suppress) {
      return;
    }

    SetErrorHandler(m_previousRootErrorHandler);
    previousRootErrorHandler = nullptr;
    genfit::errorOut.rdbuf(m_genfitErrorBuffer);
    genfit::debugOut.rdbuf(m_genfitDebugBuffer);
  }

  ScopedFitDiagnostics(const ScopedFitDiagnostics&) = delete;
  ScopedFitDiagnostics& operator=(const ScopedFitDiagnostics&) = delete;

private:
  bool m_suppress;
  NullStreamBuffer m_sink;
  std::streambuf* m_genfitErrorBuffer{nullptr};
  std::streambuf* m_genfitDebugBuffer{nullptr};
  ErrorHandlerFunc_t m_previousRootErrorHandler{nullptr};
};

} // namespace

#endif // GENFIT_DEBUG_LEVEL_H
