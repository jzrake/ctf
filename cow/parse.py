
import re
import sys

cowh = open('cow.h', 'r')
prototype = re.compile(r"\b([a-z\s]+)\s+(cow_\w+)(.*)")
macrodef = re.compile(r"#define COW_(\w+)\s+(.+)\s*")

func_proto = """
static int _%(funcname)s(lua_State *L)
{
  %(getargs)s
  %(call)s
  %(push)s
}"""

fbodies = [ ]
wrapped = [ ]
macros = [ ]

for line in cowh:
    if line.startswith('typedef') or line.startswith('struct'): continue
    if line.startswith('char *'):
        line = line.replace('char *', 'charstar ')

    if '//' in line:
        line = line[:line.index('//')-1]

    m = prototype.match(line[:-2])
    if not m:
        if line.startswith('#define'):
            n = macrodef.match(line)
            if n:
                macros.append(n.groups())
        continue
    retval, funcname, argstr = m.groups()

    args = argstr[1:-1].split(',')
    if not args[0].startswith('cow'):
        continue
    getargs = [ ]
    argnames = [ ]
    narg = 1
    for arg in args:
        argtype, argname = arg.split()[:-1], arg.split()[-1]
        argtype = ''.join(argtype)
        if argname.startswith('**'):
            argname = argname[2:]
            argtype += ' **'
        elif argname.startswith('*'):
            argname = argname[1:]
            argtype += ' *'

        ga = None

        if argtype == "cow_domain *":
            ga = """cow_domain *%s = *((cow_domain**) """\
                """luaL_checkudata(L, %d, "cow::domain"));""" % (argname, narg)
        elif argtype == "cow_dfield *":
            ga = """cow_dfield *%s = *((cow_dfield**) """\
                """luaL_checkudata(L, %d, "cow::dfield"));""" % (argname, narg)
        elif argtype == "cow_histogram *":
            ga = """cow_histogram *%s = *((cow_histogram**) """\
                """luaL_checkudata(L, %d, "cow::histogram"));""" % (argname, narg)
        elif argtype == "int":
            ga = """int %s = luaL_checkinteger(L, %d);""" % (argname, narg)
        elif argtype == "double":
            ga = """double %s = luaL_checknumber(L, %d);""" % (argname, narg)
        elif argtype == "char *":
            ga = """char *%s = (char*)luaL_checkstring(L, %d);""" % (argname, narg)
        elif argtype == "int *":
            ga = """int *%s = (int*) lua_touserdata(L, %d);""" % (argname, narg)
        elif argtype == "double *":
            ga = """double *%s = (double*) lua_touserdata(L, %d);""" % (argname, narg)
        elif argtype == "double **":
            ga = """double **%s = (double**) lua_touserdata(L, %d);""" % (argname, narg)
        elif argtype == "void *":
            ga = """void *%s = lua_touserdata(L, %d);""" % (argname, narg)
        elif argtype == "cow_transform":
            ga = """cow_transform %s = *((cow_transform*) """\
                """luaL_checkudata(L, %d, "cow::transform"));""" % (argname, narg)
        else:
            ga = """void *%s = lua_touserdata(L, %d);""" % (argname, narg)

        if ga == None: ga = "unknown " + argtype
        getargs.append(ga)
        argnames.append(argname)
        narg += 1

    call = funcname + '(' + ', '.join(argnames) + ');'

    if retval == 'void':
        push = "return 0;"
    elif retval == 'charstar':
        call = "char *ret = " + call
        push = "lua_pushstring(L, ret);\n  return 1;"
    else:
        call = ("%s ret = " % retval) + call
        push = "lua_pushnumber(L, ret);\n  return 1;"
        
    fbodies.append(func_proto % {'funcname': funcname,
                                 'getargs': '\n  '.join(getargs),
                                 'call': call,
                                 'push': push})
    wrapped.append(funcname)


wrap = open('cowfuncs.c', 'w')

for fbody in fbodies:
    wrap.write(fbody)

wrap.write("\nstatic luaL_Reg cow_module_funcs[] = {\n")
for f in wrapped:
    wrap.write("""    {"%s", _%s},\n""" % (f[4:], f))
wrap.write("    {NULL, NULL}};\n")

wrap.write("static void register_constants(lua_State *L)\n{\n")
for m in macros:
    wrap.write("""  lua_pushnumber(L, %s); lua_setfield(L, -2, "%s");\n""" % (
            m[1], m[0]))
wrap.write("}\n")
