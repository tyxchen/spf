# Install script for directory: /Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/spf" TYPE FILE FILES
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/compact_particle_population.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/csmc.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/numerical_utils.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/param.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/particle.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/particle_genealogy.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/particle_population.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/permutation_stream.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/pg.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/pg_proposal.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/pmcmc_options.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/pmmh.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/pmmh_proposal.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/resampling.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/sampling_utils.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/smc.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/smc_model.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/smc_options.hpp"
    "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/spf.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE STATIC_LIBRARY FILES "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/Debug/libSPF.a")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
    endif()
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE STATIC_LIBRARY FILES "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/Release/libSPF.a")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
    endif()
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE STATIC_LIBRARY FILES "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/MinSizeRel/libSPF.a")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
    endif()
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE STATIC_LIBRARY FILES "/Users/seonghwanjun/Dropbox/Research/smc-research/repos/spf/src/RelWithDebInfo/libSPF.a")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libSPF.a")
    endif()
  endif()
endif()

