--------------------------------------------------------------------------------
--
--        Partial implementation of the Python class model for Lua
--
--                   Copyright (c) 2012, Jonathan Zrake
--
--------------------------------------------------------------------------------
--
-- Exports the following functions:
-- 
-- + class
-- + super
-- + isclass
-- + isinstance
-- + issubclass
-- + getattrib
-- + setattrib
-- + classname
-- + mro (method resolution order)
-- + getattrib
-- + setattrib
-- 
-- and the following classes:
--
-- + object
--
-- Classes support multiple inheritance and casting via super to any base
-- class. The method resolution order is depth-first. All classes inherit from
-- object.
--------------------------------------------------------------------------------

local class_meta = { }
local instance_meta = { }
local class_module = { }

--------------------------------------------------------------------------------
-- Private utility functions
--------------------------------------------------------------------------------
local function deque()
   --
   -- Return an object which acts as a minimal double-ended queue. That is, it
   -- may be used as a queue or a stack by appropriate calls to the
   -- [push|pop]_[front|back] functions. This function is only here to aid the
   -- building of MRO's.
   --
   local self = { }
   self._items = { }
   self._front = 0
   self._back = 0
   function self:push_front(item)
      if item == nil then
	 error('cannot push nil onto the deque')
      end
      self._front = self._front - 1
      self._items[self._front] = item
   end
   function self:pop_front()
      local item = self._items[self._front]
      self._items[self._front] = nil
      if item then
	 self._front = self._front + 1
      end
      return item
   end
   function self:push_back(item)
      if item == nil then
	 error('cannot push nil onto the deque')
      end
      self._items[self._back] = item
      self._back = self._back + 1
   end
   function self:pop_back()
      local item = self._items[self._back-1]
      self._items[self._back-1] = nil
      if item then
	 self._back = self._back - 1
      end
      return item
   end
   function self:size()
      return self._back - self._front
   end
   function self:items()
      return pairs(self._items)
   end
   function self:empty()
      return not (next(self._items) and true or false)
   end
   return self
end


local function resolve(obj, key)
   -- 
   -- Attempt to resolve attribute `key` of `obj` which may be a class or
   -- instance. Recurse through __base__ table breadth-first using the MRO, and
   -- return nil if attribute is not found. Do not invoke the __index__
   -- meta-method in doing so.
   -- 
   local val = obj.__dict__[key] or obj.__class__.__dict__[key]
   if val then return val end
   for _,b in ipairs(obj.__class__.__mro__) do
      val = b.__dict__[key]
      if val then return val end
   end
   return nil
end


--------------------------------------------------------------------------------
-- Module functions
--------------------------------------------------------------------------------
local function class(name, ...)
   local base = {...}
   if #base == 0 and name ~= 'object' then
      base[1] = class_module.object
   end

   local mrovec = { }
   local mroset = { }
   local Q = deque()
   for i,b in ipairs(base) do Q:push_back(b) end

   while not Q:empty() do
      local c = Q:pop_front()
      for _,b in ipairs(c.__base__) do
	 Q:push_back(b)
      end
      if not mroset[c] then
	 mroset[c] = c
	 table.insert(mrovec, c)
      end
   end
   return setmetatable({__name__=name,
                        __base__=base,
			__mro__=mrovec,
                        __dict__={ },
			__class__={__dict__={ }, __mro__={ }}}, class_meta)
end
local function super(instance, base)
   if not base then
      return instance.__base__[1]()
   else
      for i,v in ipairs(instance.__base__) do
         if v == base then
            local proxy = v()
            rawset(proxy, '__dict__', instance.__dict__)
            return proxy
         end
      end
   end
   -- returns nil if super call cannot be made
end
local function isclass(c)
   return getmetatable(c) == class_meta
end
local function isinstance(instance, class)
   return instance.__class__ == class
end
local function issubclass(instance, base)
   return super(instance, base) and true or false
end
local function getattrib(instance, key)
   return resolve(instance, key)
end
local function setattrib(instance, key, value)
   instance.__dict__[key] = value
end
local function classname(A)
   return isclass(A) and A.__name__ or A.__class__.__name__
end
local function mro(A)
   return A.__mro__
end

--------------------------------------------------------------------------------
-- Class metatable
--------------------------------------------------------------------------------
function class_meta:__call(...)
   local dict = { }
   local base = { }
   for k,v in pairs(self.__base__) do base[k] = v end
   local new = setmetatable({__name__=self.__name__,
                             __dict__=dict,
                             __class__=self,
			     __base__=self.__base__}, instance_meta)
   if getattrib(new, '__init__') then
      getattrib(new, '__init__')(new, ...)
   end
   return new
end
function class_meta:__index(key)
   return resolve(self, key)
end
function class_meta:__newindex(key, value)
   self.__dict__[key] = value
end
function class_meta:__tostring(key, value)
   return string.format('<Class: %s @ %s>', classname(self),
                        string.sub(tostring(self.__dict__), 8))
end

--------------------------------------------------------------------------------
-- Instance metatable
--------------------------------------------------------------------------------
function instance_meta:__index(key)
   local index = resolve(self, '__index__')
   local def = resolve(self, key)
   if type(index) == 'function' then return def or index(self, key)
   else return def or index[key]
   end
end
function instance_meta:__newindex(key, value)
   if self.__dict__[key] then
      self.__dict__[key] = value
   else
      resolve(self, '__newindex__')(self, key, value)
   end
end
function instance_meta:__tostring()
   return resolve(self, '__tostring__')(self)
end
function instance_meta:__pairs()
   return resolve(self, '__pairs__')(self)
end
function instance_meta:__gc()
   return resolve(self, '__gc__')(self)
end

--------------------------------------------------------------------------------
-- object Class definition
--------------------------------------------------------------------------------
local object = class('object')
function object:__index__(key)
   return resolve(self, key)
end
function object:__newindex__(key, value)
   self.__dict__[key] = value
end
function object:__tostring__()
   return string.format('<Class instance: %s @ %s>', classname(self),
                        string.sub(tostring(self.__dict__), 8))
end
function object:__pairs__()
   error('object does not support iteration')
end
function object:__call__()
   error('object does not support call')
end
function object:__gc__()
   -- warning! __gc__ is only triggered for the first resolved __gc__ method in
   -- the hierarchy.
end

--------------------------------------------------------------------------------
-- Class module definition
--------------------------------------------------------------------------------
class_module.class = class
class_module.super = super
class_module.isclass = isclass
class_module.isinstance = isinstance
class_module.issubclass = issubclass
class_module.getattrib = getattrib
class_module.setattrib = setattrib
class_module.classname = classname
class_module.mro = mro
class_module.object = object

--------------------------------------------------------------------------------
-- Unit test
--------------------------------------------------------------------------------
local function test()
   local SoftObject = class('SoftObject')
   function SoftObject:set_softness(val)
      self._softness = val
   end
   function SoftObject:get_softness(val)
      return self._softness
   end

   local Animal = class('Animal')
   function Animal:speak()
      return 'unknown noise'
   end
   function Animal:eat()
      return 'unknown food'
   end
   function Animal:jump()
      return 'cannot jump'
   end

   local Cat = class('Cat', Animal, SoftObject)
   function Cat:__init__(softness)
      self._softness = softness
   end
   function Cat:__tostring__()
      return '<:crazy cat:>'
   end
   function Cat:__index__(key)
      if key == 'food' then
         return 'starving'
      end
   end
   function Cat:__newindex__(key, value) -- boring over-ride
      self.__dict__[key] = value
   end
   function Cat:speak()
      return 'meow'
   end
   function Cat:jump()
      return 'can jump'
   end
   local blue = Cat(100)

   blue.tree = 'blue tree'
   assert(blue.food == 'starving')
   assert(blue.tree == 'blue tree')
   assert(blue:jump() == 'can jump')
   assert(type(blue.get_softness) == 'function')
   assert(blue:get_softness() == 100)
   assert(blue:speak() == 'meow')
   assert(blue:eat() == 'unknown food')
   assert(super(blue):speak() == 'unknown noise')
   assert(super(blue, Animal):jump() == 'cannot jump')
   assert(tostring(blue) == '<:crazy cat:>')

   -- proxy class returned by super retains __dict__
   assert(super(blue, Animal).tree == 'blue tree')

   assert(isclass(Animal))
   assert(isinstance(blue, Cat))
   assert(issubclass(blue, Animal))
   assert(issubclass(blue, SoftObject))
   assert(not isclass({ }))
   assert(not isclass(blue))
   assert(not issubclass(blue, class('BogusBase')))
end

--------------------------------------------------------------------------------
-- Run test or export module
--------------------------------------------------------------------------------
if ... then -- if __name__ == "__main__"
   return class_module
else
   test()
   print(debug.getinfo(1).source, ": All tests passed")
end
