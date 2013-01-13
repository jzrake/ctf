
#include "logging.hpp"
#include <iostream>

DebugLogger::DebugLevel DebugLogger::Level = DebugLogger::LevelInfo;
DebugLogger DebugLog;

void DebugLogger::erase()
{
  LogBuffer.str(std::string());
  NullStream.str(std::string());
}
std::string DebugLogger::str() const
{
  return LogBuffer.str();
}

std::ostream &DebugLogger::Info(std::string caller)
{
  if (Level <= LevelInfo) {
    LogBuffer << "Info from " << caller << ": " << std::endl;
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Info()
{
  if (Level <= LevelInfo) {
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Warning(std::string caller)
{
  if (Level <= LevelWarning) {
    LogBuffer << "Warning from " << caller << ": " << std::endl;
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Warning()
{
  if (Level <= LevelWarning) {
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Error(std::string caller)
{
  if (Level <= LevelError) {
    LogBuffer << "Error from " << caller << ": " << std::endl;
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Error()
{
  if (Level <= LevelError) {
    return LogBuffer;
  }
  else
    return NullStream;
}
std::ostream &DebugLogger::Panic(std::string caller)
{
  std::cout << "Panic from " << caller << ": " << std::endl;
  return std::cout;
}
std::ostream &DebugLogger::Panic() {
  return std::cout;
}
