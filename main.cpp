//#include <iostream>
#include <cstdint>
#include <iostream>
#include <print>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <raylib.h>


uint64_t window_width = 900;
uint64_t window_height = 600;

double target_fps = 60.f;

double aspect_ratio = 16.f / 9.f;

uint64_t image_width = 400;
uint64_t image_height = image_width / aspect_ratio;

double viewport_height = 2.f;
double viewport_width = viewport_height * (double)image_width / image_height;

double focal_length = 1.f;


struct Vec3 {
    double x;
    double y;
    double z;
	
	double length() const {
	    return sqrtf(x*x + y*y + z*z);
	}
};

struct Sphere {
    Vec3 center;
    double r;
};

Sphere sphere = {{0, 0, -1}, 0.5f};

void add_inplace(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

Vec3 operator+(const Vec3& a, const Vec3& b) {
    Vec3 v;
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    Vec3 v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;
    return v;
}

Vec3 operator*(const Vec3& v, double t) {
    Vec3 res;
    res.x = v.x * t;
    res.y = v.y * t;
    res.z = v.z * t;
    return res;
}

Vec3 operator/(const Vec3& v, double t) {
    Vec3 res;
    res.x = v.x / t;
    res.y = v.y / t;
    res.z = v.z / t;
    return res;
}

std::ostream& operator<<(std::ostream& outstr, const Vec3 v) {
	outstr << v.x << " " << v.y << " " << v.z; 
	return outstr;
}

std::string operator<<(const Vec3 v, const char* str) {
	return str;
}

Vec3 normalize(const Vec3& v) {
	return v / v.length();
}

double dot(Vec3 a, Vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 lerp(Vec3 start, Vec3 end, double t) {
	return start * (1.f - t) + end * t; 
}

Color lerp(Color start, Color end, double t) {
	double inv_t = 1.f - t;
	Color res = {
			(uint8_t)(start.r * inv_t + end.r * t), 
			(uint8_t)(start.g * inv_t + end.g * t),
			(uint8_t)(start.b * inv_t + end.b * t),
			255
	};
	return res;
}

struct RayD {
    Vec3 origin;
    Vec3 dir;

    Vec3 at(double t) const {
	return origin + dir * t;
    }
};

double hit_sphere(const Sphere& sphere, const RayD& ray) {
	Vec3 C_Q = sphere.center - ray.origin;
	double a = dot(ray.dir, ray.dir);
	//double b = dot((ray.dir * -2.f), C_Q);
	double h = dot(ray.dir, C_Q);
	double c = dot(C_Q, C_Q) - sphere.r * sphere.r;

	double discriminant = h*h - a*c;

	if (discriminant < 0.f) return -1.f;
	else {
	    return (h - std::sqrtf(discriminant)) / a;
	}
}

Color ray_color(const RayD& ray) {
	double t = hit_sphere(sphere, ray);

	if (t > 0) {

		Vec3 hit = ray.at(t);
		Vec3 norm_s = normalize(hit - sphere.center);

		return {.r = (uint8_t)(255.f * (norm_s.x + 1.f) / 2.f), 
			.g = (uint8_t)(255.f * (norm_s.y + 1.f) / 2.f), 
			.b = (uint8_t)(255.f * (norm_s.z + 1.f) / 2.f),
			.a = 255};
	}
	else {

		// background
		Color end = {127, (uint8_t)(0.7f * 255), 255};
		Color start = {255, 255, 255};

		Vec3 norm = normalize(ray.dir);
		// y could be -1 < y < 1 so we bring it to 0 - 1
		double t2 = (norm.y + 1.f) / 2.f;
		assert(t2 <= 1.f && t2 >= 0.f && "t is not in the form 0 <= t <= 1");
		return lerp(start, end, t2);
	}
}

// assumes a is [0..1+]
uint8_t float_to_byte(double a) {
    return (uint8_t) (a * 255.999f);
}

int main() {

	InitWindow(window_width, window_height, "Raytracer");
	SetTargetFPS(target_fps);

	Image img = GenImageColor(image_width, image_height, BLACK);

	Vec3 camera_pos = {0};


	Vec3 pixel_delta_u = {viewport_width / image_width, 0, 0};
	Vec3 pixel_delta_v = {0, -viewport_height / image_height, 0};
	
	Vec3 viewport_top_left = {-viewport_width / 2.f, viewport_height / 2.f, -focal_length}; 
	Vec3 pixel_start_pos = viewport_top_left + pixel_delta_u / 2.f + pixel_delta_v / 2.f;

	Texture txt = LoadTextureFromImage(img);

	while (!WindowShouldClose()) {
		float step = 1.f / target_fps;
		if (IsKeyDown(KEY_LEFT)) {
			camera_pos.x -= step;
		}
		if (IsKeyDown(KEY_RIGHT)) {
			camera_pos.x += step;
		}
		if (IsKeyDown(KEY_UP)) {
			camera_pos.z -= step;
		}
		if (IsKeyDown(KEY_DOWN)) {
			camera_pos.z += step;
		}

		BeginDrawing();

		for (int y = 0; y < image_height; ++y) {
		    for (int x = 0; x < image_width; ++x) {

			Vec3 pixel_center = pixel_start_pos + pixel_delta_u * x + pixel_delta_v * y;
			Vec3 ray_dir = pixel_center - camera_pos;
			
			RayD ray = {camera_pos, ray_dir};
			
			Color c = ray_color(ray);
			ImageDrawPixel(&img, x, y, c);
		    }
		}

		UpdateTexture(txt, img.data);
		DrawTexture(txt, 0, 0, WHITE);
		EndDrawing();
	}

	CloseWindow();
	
	//std::println("P3");
	//std::println("{} {}", image_width, image_height);
	//std::println("255");
	//
	//Vec3 camera_pos = {0};


	//Vec3 pixel_delta_u = {viewport_width / image_width, 0, 0};
	//Vec3 pixel_delta_v = {0, -viewport_height / image_height, 0};
	//
	//Vec3 viewport_top_left = {-viewport_width / 2.f, viewport_height / 2.f, -focal_length}; 
	//Vec3 pixel_start_pos = viewport_top_left + pixel_delta_u / 2.f + pixel_delta_v / 2.f;

	//for (int y = 0; y < image_height; ++y) {
	//    std::clog << "\rScanlines remaining: " << (image_height - y) << std::flush;
	//    //Color c(0, float_to_byte((double)y / (image_height - 1)), 0.f);
	//    for (int x = 0; x < image_width; ++x) {
	//	//c.r = float_to_byte((double)x / (image_width - 1));

	//	Vec3 pixel_center = pixel_start_pos + pixel_delta_u * x + pixel_delta_v * y;
	//	Vec3 ray_dir = pixel_center - camera_pos;
	//	
	//	RayD ray = {camera_pos, ray_dir};
	//	
	//	Color c = ray_color(ray);
	//	std::println("{} {} {}", c.r, c.g, c.b);
	//    }
	//}
	//std::clog << "\nDone";
	

	return 0;
}

