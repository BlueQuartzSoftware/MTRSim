/**
 * The LibMTRSim sources are compiled directly into the SIMPLNX MTRSim plugin
 * shared library (rather than a standalone libMTRSim). This shim defines the
 * LibMTRSim export macro so the LibMTRSim symbols are exported from / imported
 * into the plugin DLL exactly like the plugin's own symbols. Without this, on
 * MSVC the LibMTRSim symbols (e.g. mtrsim::readODFComponents) are not exported
 * from the plugin DLL and external consumers -- the unit-test executable and
 * the Python bindings -- fail to link with LNK2019 / LNK1120. On non-Windows
 * the symbols use default visibility (the prior behavior).
 *
 * CMake defines `MTRSim_EXPORTS` while building the plugin target.
 */
#pragma once

#ifndef LIBMTRSIM_EXPORT
#ifdef _WIN32
#ifdef MTRSim_EXPORTS
#define LIBMTRSIM_EXPORT __declspec(dllexport)
#else
#define LIBMTRSIM_EXPORT __declspec(dllimport)
#endif
#else
#define LIBMTRSIM_EXPORT __attribute__((visibility("default")))
#endif
#endif
