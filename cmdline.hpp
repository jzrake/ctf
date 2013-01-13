

#ifndef __CommandLineParser_HEADER__
#define __CommandLineParser_HEADER__

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <list>

namespace Opt
{
  enum OptionType { Bool, Integer, Double, String };
}

class OptionParser
{
public:
private:
  struct CommandLineOption
  {
    CommandLineOption(void *addr, Opt::OptionType type, char snam,
                      const char *lnam, const char *help)
      : addr(addr), type(type), snam(snam), lnam(lnam), help(help), used(0) { }

    void *addr;
    Opt::OptionType type;
    const char snam, *lnam, *help;
    int used;
  } ;
  std::list<CommandLineOption> RegisteredOptions;

public:
  class DuplicateOption : public std::exception
  {
  private:
    const char *Key;
  public:
    DuplicateOption(const char *Key) : Key(Key) { }
    virtual const char *what() const throw()
    {
      return Key;
    }
  } ;
  class UnrecognizedOption : public std::exception
  {
  private:
    const char *Key;
  public:
    UnrecognizedOption(const char *Key) : Key(Key) { }
    virtual const char *what() const throw()
    {
      return Key;
    }
  } ;
  class OptionRequiresValue : public std::exception
  {
  private:
    const char *Key;
  public:
    OptionRequiresValue(const char *Key) : Key(Key) { }
    virtual const char *what() const throw()
    {
      return Key;
    }
  } ;
  class OptionForbidsValue : public std::exception
  {
  private:
    const char *Key;
  public:
    OptionForbidsValue(const char *Key) : Key(Key) { }
    virtual const char *what() const throw()
    {
      return Key;
    }
  } ;

  int HelpMessage;
  OptionParser() : HelpMessage(0)
  {
    Register(&HelpMessage, Opt::Bool, 'h', "help", "Display help message");
  }

  void Register(void *addr, Opt::OptionType type, char snam,
                const char *lnam,
                const char *help)
  {
    std::list<CommandLineOption>::const_iterator it = RegisteredOptions.begin();
    while (it != RegisteredOptions.end()) {
      if (it->snam == snam) throw DuplicateOption(&snam);
      if (!strcmp(it->lnam, lnam)) throw DuplicateOption(lnam);
      ++it;
    }
    CommandLineOption opt(addr, type, snam, lnam, help);
    RegisteredOptions.push_back(opt);
  }

  void Register(void *addr, Opt::OptionType type,
                const char *lnam,
                const char *help)
  {
    std::list<CommandLineOption>::const_iterator it = RegisteredOptions.begin();
    while (it != RegisteredOptions.end()) {
      if (!strcmp(it->lnam, lnam)) throw DuplicateOption(lnam);
      ++it;
    }
    CommandLineOption opt(addr, type, '?', lnam, help);
    RegisteredOptions.push_back(opt);
  }
  void PrintUsageMessage()
  {
    std::list<CommandLineOption>::const_iterator it = RegisteredOptions.begin();
    while (it != RegisteredOptions.end()) {
      if (it->snam == '?')
        printf("--%s: %s\n", it->lnam, it->help);
      else
        printf("--%s -%c: %s\n", it->lnam, it->snam, it->help);
      ++it;
    }
  }
  std::vector<std::string> ParseOptions(int argc, char **argv)
  {
    std::vector<std::string> args;

    for (int i=1; i<argc; ++i) {
      const char *ai = argv[i];

      if      (ai[0] == '-' && ai[1] != '-') {
        ProcessShortOption(ai[1], ai+2);
      }
      else if (ai[0] == '-' && ai[1] == '-') {
        i += ProcessLongOption(ai+2, argv[i+1]);
      }
      else {
        args.push_back(ai);
      }
    }
    return args;
  }
private:
  void ProcessShortOption(char snam, const char *val)
  {
    std::list<CommandLineOption>::iterator it = RegisteredOptions.begin();
    while (it != RegisteredOptions.end()) {
      if (it->snam == snam) {
        switch (it->type) {

        case Opt::Bool:
          if (it->used) throw DuplicateOption(&snam);
          if (*val) throw OptionForbidsValue(val);
          *(int*)it->addr = 1;
          it->used = true;
          return;

        case Opt::Integer:
          if (it->used) throw DuplicateOption(&snam);
          if (!*val) throw OptionRequiresValue(&snam);
          *(int*)it->addr = atoi(val);
          it->used = true;
          return;

        case Opt::Double:
          if (it->used) throw DuplicateOption(&snam);
          if (!*val) throw OptionRequiresValue(&snam);
          *(double*)it->addr = atof(val);
          it->used = true;
          return;

        case Opt::String:
          if (it->used) throw DuplicateOption(&snam);
          if (!*val) throw OptionRequiresValue(&snam);
          *(const char**)it->addr = val;
          it->used = true;
          return;
        }
      }
      ++it;
    }
    throw UnrecognizedOption(&snam);
  }
  int ProcessLongOption(const char *lnam, const char *val)
  {
    std::list<CommandLineOption>::iterator it = RegisteredOptions.begin();

    while (it != RegisteredOptions.end()) {
      if (!strcmp(it->lnam, lnam)) {
        switch (it->type) {

        case Opt::Bool:
          if (it->used) throw DuplicateOption(lnam);
          *(int*)it->addr = 1;
          it->used = true;
          return 0;

        case Opt::Integer:
          if (it->used) throw DuplicateOption(lnam);
          if (val[0] == '-') throw OptionRequiresValue(lnam);
          *(int*)it->addr = atoi(val);
          it->used = true;
          return 1;

        case Opt::Double:
          if (it->used) throw DuplicateOption(lnam);
          if (val[0] == '-') throw OptionRequiresValue(lnam);
          *(double*)it->addr = atof(val);
          it->used = true;
          return 1;

        case Opt::String:
          if (it->used) throw DuplicateOption(lnam);
          if (val[0] == '-') throw OptionRequiresValue(lnam);
          *(const char**)it->addr = val;
          it->used = true;
          return 1;
        }
      }
      ++it;
    }
    throw UnrecognizedOption(lnam);
  }
} ;


#endif // __CommandLineParser_HEADER__
