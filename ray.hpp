#pragma once
#include <gmath.hpp>

struct RayD {
    Vec3 origin;
    Vec3 dir;

    Vec3 at(double t) const {
	return origin + dir * t;
    }
};
