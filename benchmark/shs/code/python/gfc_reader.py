"""Minimal ICGEM .gfc header reader used to validate the shared model before any engine reads it.

Scope (deliberately small):
  * This module reads ONLY the header block (up to ``end_of_head``) and checks the constants
    the benchmark depends on: GM, reference radius, max degree, normalisation and tide system.
  * It does NOT parse the coefficients. Each engine reads the .gfc with its own reader
    (pyshtools: ``SHGravCoeffs.from_file(format='icgem')``; pyharm: ``Shc.from_file('gfc')``).
    The runners then compare the GM and R that each library reports against this header.
  * The tide system is checked here because the pyshtools ICGEM reader ignores it
    (its docstring says ``tide_system (not used)``). The benchmark requires ``tide_free``.

Header keys follow the ICGEM GFC specification (``earth_gravity_constant``, ``radius``,
``max_degree``, ``norm``, ``tide_system``). Fortran exponents (``D``) are accepted.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class GfcHeader:
    """Constants read from the header of an ICGEM .gfc file (SI units)."""

    path: str
    modelname: str
    gm: float          # gravitational constant times mass of the Earth [m^3/s^2]
    radius: float      # reference radius of the coefficients [m]
    max_degree: int    # maximum degree stored in the file
    norm: str          # e.g. 'fully_normalized'
    tide_system: str   # must be 'tide_free' for this benchmark


def _to_float(text):
    return float(text.strip().replace('D', 'E').replace('d', 'e'))


def read_gfc_header(path):
    """Read the header of a .gfc file and return a :class:`GfcHeader`.

    Raises ValueError if a required key is missing, if the normalisation is not fully
    normalised, or if the file is not tide-free.
    """
    keys = {}
    with open(path, 'r', encoding='ascii', errors='replace') as fh:
        for line in fh:
            text = line.strip()
            if not text:
                continue
            if text.startswith('end_of_head'):
                break
            parts = text.split(None, 1)
            if len(parts) == 2:
                keys[parts[0]] = parts[1].strip()
        else:
            raise ValueError(f'{path}: no end_of_head line found')

    def need(name):
        if name not in keys:
            raise ValueError(f'{path}: header key {name!r} is missing')
        return keys[name]

    gm_key = 'earth_gravity_constant' if 'earth_gravity_constant' in keys else 'gravity_constant'
    header = GfcHeader(
        path=str(path),
        modelname=keys.get('modelname', 'unknown'),
        gm=_to_float(need(gm_key)),
        radius=_to_float(need('radius')),
        max_degree=int(need('max_degree')),
        norm=need('norm').lower(),
        tide_system=need('tide_system').lower(),
    )
    if header.norm != 'fully_normalized':
        raise ValueError(f'{path}: norm is {header.norm!r}, expected fully_normalized')
    if header.tide_system != 'tide_free':
        raise ValueError(f'{path}: tide_system is {header.tide_system!r}, expected tide_free')
    return header
