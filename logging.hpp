
#ifndef __DebugLogger_HEADER__
#define __DebugLogger_HEADER__

#include <sstream>

class DebugLogger
{
public:
  enum DebugLevel { LevelInfo, LevelWarning, LevelError, LevelPanic };
  static DebugLevel Level;

private:
  std::stringstream LogBuffer, NullStream;

public:
  DebugLogger() { }
  ~DebugLogger() { }

  void erase();
  std::string str() const;

  std::ostream &Info(std::string caller);
  std::ostream &Info();
  std::ostream &Warning(std::string caller);
  std::ostream &Warning();
  std::ostream &Error(std::string caller);
  std::ostream &Error();
  std::ostream &Panic(std::string caller);
  std::ostream &Panic();
} ;

extern DebugLogger DebugLog;

#endif // __DebugLogger_HEADER__
