#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <unistd.h>

using namespace std;

/** global variables**/
double kEps = 0.001;
double kd = 1.0;
double ls = 40.0;
double invPI = 1/M_PI;
int max_depth = 10;
bool separate = true;
int sample_per_pixel = 100;
int area_light_num_rays = 100;
int M = 256;
int N = 64;
/**
 * Vec3 object
 **/
struct Vec3 {
	double x;
	double y;
	double z;

	Vec3() {
		x = 0;
		y = 0;
		z = 0;
	}

	Vec3(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
	}

	Vec3 operator *(const double t) const {
		return Vec3(x*t, y*t, z*t);
	}

	Vec3 operator /(const double t) const {
		return Vec3(x/t, y/t, z/t);
	}

	Vec3 operator +(const Vec3 & other) const {
		return Vec3(x+other.x, y+other.y, z+other.z);
	}

	Vec3 operator -(const Vec3 & other) const {
		return Vec3(x - other.x, y-other.y, z-other.z);
	}

	double dot(const Vec3 & other) const {
		return x*other.x + y*other.y +z*other.z;
	}

	Vec3 cross(const Vec3 & other) const {
		return Vec3(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);
	}

	double len() const {
		return sqrt(x*x + y*y + z*z);
	}
};

ostream& operator <<(ostream& os, const Vec3 v) {
	os << v.x << " " << v.y << " " << v.z;
	return os;
}

/** Ray object **/
struct Ray {
	Vec3 o;
	Vec3 d;

	Ray() {}

	Ray(Vec3 a, Vec3 b) {
		o = a;
		d = b;
	}

	Vec3 getPt(double t) {
		return o+d*t;
	}
};

/** Rectangle object **/
struct Rectangle {
	Vec3 p0;
	Vec3 a;
	Vec3 b;
	Vec3 n;
	double als;
	double bls;
	Vec3 material;

	Rectangle(Vec3 p0, Vec3 a, Vec3 b) {
		this -> p0 = p0;
		this -> a = a;
		this -> b = b;
		n = b.cross(a);
		n = n / n.len();
		material = Vec3(255,0,0);
		als = a.dot(a);
		bls = b.dot(b);
	}

	Rectangle(Vec3 p0, Vec3 a, Vec3 b, Vec3 m) {
		this -> p0 = p0;
		this -> a = a;
		this -> b = b;
		n = a.cross(b);
		n = n / n.len();
		als = a.dot(a);
		bls = b.dot(b);
		material = m;
		//cout << n << endl;
	}

	bool hit(Ray & ray, double & tmin) const {
		float t = (p0 - ray.o).dot(n) / (ray.d.dot(n));
		if (t <= kEps)
			return false;
		Vec3 p = ray.o + ray.d * t;
		Vec3 d = p - p0;
		double ddota = d.dot(a);
		if (ddota < 0.0 || ddota > als)
			return false;
		double ddotb = d.dot(b);
		if (ddotb < 0.0 || ddotb > bls)
			return false;

		tmin = t;
		return true;
	}
};

/** Viewport object **/
struct Viewport {
	int w;
	int h;
	double s;
	Viewport(int a, int b, double c) {
		w = a;
		h = b;
		s = c;
	}

	Vec3 getPixelCenter(int c, int r, int d) {
		return Vec3(s * (c - w/2 + 0.5), s * (r - h/2 + 0.5),d);
	}
};

struct Irradiance {
	Vec3 p;
	Vec3 n;
	Vec3 E;
	float R;
	Irradiance(Vec3 a, Vec3 b, Vec3 c) {
		p = a;
		n = b;
		E = c;
	}
};

vector<Rectangle> recs;
vector<double> samples;
vector<Irradiance> cache;

Vec3 color_multiply(Vec3 c1, Vec3 c2) {
	Vec3 res = Vec3(0,0,0);
	res.x = int(c1.x * c2.x / 255);
	res.y = int(c1.y * c2.y / 255);
	res.z = int(c1.z * c2.z / 255);
	return res;
}

Vec3 clamp(Vec3 & color) {
	color.x = max(0,min(255,(int)(color.x)));
	color.y = max(0,min(255,(int)(color.y)));
	color.z = max(0,min(255,(int)(color.z)));
	return color;
}

void multi_jitter(int n) {
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			samples.push_back((i + (j + drand48()) / n) / n);
			samples.push_back((j + (i + drand48()) / n) / n);
		}
	}

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			int k = j + drand48() *(n - j);
			swap(samples[(j*n+i)*2], samples[(k*n+i)*2]);
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int k = i + drand48() * (n -i);
			swap(samples[(j*n+i)*2+1], samples[(j*n+k)*2+1]);
		}
	}
	return;
}

Vec3 area_light_shade(Vec3 pt, int idx, int num_rays) {
	Vec3 color = Vec3(0,0,0);
	Ray shadow = Ray();
	double rec_area = recs[0].a.len() * recs[0].b.len();
	//cout << rec_area << endl;
	for (int i = 0; i < num_rays; i++) {
		Vec3 sample = recs[0].p0 + recs[0].a * drand48() + recs[0].b * drand48();
		Vec3 wi = sample - pt;
		double d = wi.len();
		double invd2 = 1/(d*d);
		//cout << d*d << endl;
		wi = wi / d;
		double ndotwi = recs[idx].n.dot(wi);
		if (ndotwi > 0.0) {
			bool in_shadow = false;
			shadow.o = pt;
			shadow.d = wi;
			double t_shadow;
			for (int j = 1; j < recs.size(); j++) {
				if (recs[j].hit(shadow, t_shadow) && t_shadow > 0.1 && t_shadow < d) {
					in_shadow = true;
					break;
				}
			}
			if (!in_shadow) {
				double ndotd = recs[0].n.dot(wi*(-1));
				color = color + color_multiply(recs[idx].material, recs[0].material) * kd * invPI * ls * ndotwi * ndotd * invd2 * rec_area/ num_rays;
			}
		}
	}
	clamp(color);
	//cout << color <<endl;
	return color;
}

Vec3 sample_f(const Vec3& wo, Vec3& wi, float& pdf, Vec3 pt, int idx) {
	Vec3 w = recs[idx].n;
	Vec3 v = Vec3(0.0021, 1.0, 0.0031).cross(w);
	v = v / v.len();
	Vec3 u = v.cross(w);

	double x = drand48();
	float cos_phi = cos(2 * M_PI * x);
	float sin_phi = sin(2 * M_PI * x);
	float cos_theta = 1.0 - drand48();
	float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	float pu = sin_theta * cos_phi;
	float pv = sin_theta * sin_phi;
	float pw = cos_theta;
	wi = u * pu + v * pv + w * pw;
	wi = wi / wi.len();
	//cout << wi << endl;
	pdf = recs[idx].n.dot(wi) * invPI;
	return recs[idx].material * kd * invPI;
}

Vec3 sample_f(const Vec3& wo, Vec3& wi, float& pdf, Vec3 pt, int idx, int j, int k) {
	Vec3 w = recs[idx].n;
	Vec3 v = Vec3(0.0021, 1.0, 0.0031).cross(w);
	v = v / v.len();
	Vec3 u = v.cross(w);
 	float x = drand48();
	float y = drand48();
	float cos_phi = cos(2 * M_PI / M * (j + x));
	float sin_phi = sin(2 * M_PI / M * (j + x));
	float cos_theta = 1.0 - 1.0 / N * (k + y);
	float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	float pu = sin_theta * cos_phi;
	float pv = sin_theta * sin_phi;
	float pw = cos_theta;
	wi = u * pu + v * pv + w * pw;
	wi = wi / wi.len();
	//cout << wi << endl;
	pdf = recs[idx].n.dot(wi) * invPI;
	return recs[idx].material * kd * invPI;
}

Vec3 path_shade(const Ray ray, const int depth, Vec3 pt, int idx, long long& ray_count) {
	Vec3 L = Vec3(0,0,0);
	// if (depth == 0 && separate) {
	// 	ray_count += area_light_num_rays;
	// 	L = area_light_shade(pt, idx, area_light_num_rays);
	// }
	Vec3 wi;
	Vec3 wo = ray.d * (-1);
	float pdf;
	Vec3 f = sample_f(wo, wi, pdf, pt, idx);
	float ndotwi = recs[idx].n.dot(wi);
	Ray reflected_ray(pt, wi);
	if (depth + 1 > max_depth)
		return L;
	else {
		ray_count ++;
		Vec3 next_hit;
		int next_idx;
		vector<double> t_arr;
		for (int i = 0; i < recs.size(); i++) {
			double t;
			if (recs[i].hit(reflected_ray,t) == true) {
				t_arr.push_back(t);
			}
			else {
				t_arr.push_back(DBL_MAX);
			}
		}
		vector<double>::iterator it = min_element(t_arr.begin(),t_arr.end());
		next_idx = it - t_arr.begin();
		double t = t_arr[next_idx];

		if (t != DBL_MAX) {
			//hit an object
			if (next_idx == 0) {
				//emissive
				if (depth == 0 && separate) {
					return L;
				}
				else {
					L = L + color_multiply(f, recs[0].material * 20) * ndotwi/pdf;
					return L;
				}
			}
			else {
				//matt
				next_hit = reflected_ray.getPt(t);
				L = L + color_multiply(f, path_shade(reflected_ray, depth + 1, next_hit, next_idx, ray_count)) * ndotwi /pdf;
				return L;
			}
		}
		else {
			//escape scene
			return L;
		}
	}
}

Vec3 path_shade(const Ray& ray, const int depth, Vec3 pt, int idx, long long& ray_count, int j, int k, float& inv_d, double& min_R) {
	Vec3 L = Vec3(0,0,0);
	// if (depth == 0 && separate) {
	// 	ray_count += area_light_num_rays;
	// 	L = area_light_shade(pt, idx, area_light_num_rays);
	// }
	Vec3 wi;
	Vec3 wo = ray.d * (-1);
	float pdf;
	Vec3 f = sample_f(wo, wi, pdf, pt, idx, j, k);
	float ndotwi = recs[idx].n.dot(wi);
	Ray reflected_ray(pt, wi);
	if (depth + 1 > max_depth)
		return L;
	else {
		ray_count ++;
		Vec3 next_hit;
		int next_idx;
		vector<double> t_arr;
		for (int i = 0; i < recs.size(); i++) {
			double t;
			if (recs[i].hit(reflected_ray,t) == true) {
				t_arr.push_back(t);
			}
			else {
				t_arr.push_back(DBL_MAX);
			}
		}
		vector<double>::iterator it = min_element(t_arr.begin(),t_arr.end());
		next_idx = it - t_arr.begin();
		double t = t_arr[next_idx];

		if (t != DBL_MAX) {
			inv_d = 1 / (reflected_ray.getPt(t) - pt).len();
			if ((reflected_ray.getPt(t) - pt).len() < min_R)
				min_R = (reflected_ray.getPt(t) - pt).len();
			//hit an object
			if (next_idx == 0) {
				//emissive
				if (depth == 0 && separate) {
					return L;
				}
				else {
					L = L + color_multiply(f, recs[0].material * 20) * ndotwi/pdf;
					return L;
				}
			}
			else {
				//matt
				next_hit = reflected_ray.getPt(t);
				L = L + color_multiply(f, path_shade(reflected_ray, depth + 1, next_hit, next_idx, ray_count)) * ndotwi /pdf;
				return L;
			}
		}
		else {
			inv_d = 0.0;
			//escape scene
			return L;
		}
	}
}

vector<int> bruteforceLookup(Vec3 pt, Vec3 normal) {
	vector<int> Sp;
	for (int i = 0; i < cache.size(); i++) {
		if ((pt - cache[i].p).len() < cache[i].R && normal.dot(cache[i].n) > 0.95 && (pt - cache[i].p).dot(normal + cache[i].n)/2 > -0.01) {
			if (1 / ((pt - cache[i].p).len() / cache[i].R + sqrt(1 - normal.dot(cache[i].n))) > 0.1)
				Sp.push_back(i);
		}
	}
	return Sp;
}

Vec3 irradiance_caching(const Ray& ray, Vec3 pt, Vec3 normal, int idx, long long& ray_count, bool& create, clock_t& time) {
	vector<int> Sp = bruteforceLookup(pt, normal);
	if (!Sp.empty()) {
		Vec3 color = Vec3(0,0,0);
		float w_denom = 0;
		for (int i = 0; i < Sp.size(); i++) {
			float w = 1 / ((pt - cache[Sp[i]].p).len() / cache[Sp[i]].R + sqrt(1 - normal.dot(cache[Sp[i]].n)));
			if (w > 0.1) {
				w_denom += w;
				color = color + cache[Sp[i]].E * w;
			}
		}
		color = color / w_denom;
		clamp(color);
		return color;
	}
	else {
		clock_t startTime = clock();
		create = true;
		Vec3 color = Vec3(0,0,0);
		float invD = 0;
		int count = 0;
		double min_R = DBL_MAX;
		for (int j = 0; j < M; j++)
			for (int k = 0; k < N; k++) {
				float inv_d;
				color = color +  path_shade(ray, 0, pt, idx, ray_count, j, k,inv_d, min_R);
				if (inv_d == 0.0)
					count ++;
				invD += inv_d;
			}
		color = color / (M*N);
		clamp(color);
		Irradiance record(pt,normal,color);
		//record.R = (M* N-count) / invD;
		record.R = min_R;
		// if (record.R < 1.0)
		// cout << count << endl;
		cache.push_back(record);
		time += clock() - startTime;
		return record.E;
	}
}

int main() {
	clock_t time = clock() - clock();
	long long ray_count = 0;

	Vec3 eye = Vec3(278,273,-800);
	Vec3 dir = Vec3(0,0,1);
	Vec3 up = Vec3(0,1,0);

	//Light
	Vec3 b0 = Vec3(130,0.0,0.0);
	Vec3 a0 = Vec3(0.0,0.0,-105.0);
	Vec3 p0_0 = Vec3(213.0,548.8,332.0);
	recs.push_back(Rectangle(p0_0,a0,b0, Vec3(255,186,102)));
	//left wall
	Vec3 a1 = Vec3(0.0,0.0,559.2);
	Vec3 b1 = Vec3(0.0,548.8,0.0);
	Vec3 p0_1 = Vec3(556.0,0.0,0.0);
	recs.push_back(Rectangle(p0_1,a1,b1, Vec3(255,0,0)));
	//right wall
	Vec3 b2 = Vec3(0.0,0.0,559.2);
	Vec3 a2 = Vec3(0.0,548.8,0.0);
	Vec3 p0_2 = Vec3(0.0,0.0,0.0);
	recs.push_back(Rectangle(p0_2,a2,b2, Vec3(0,255,0)));
	//back wall
	Vec3 b3 = Vec3(0,548.8,0.0);
	Vec3 a3 = Vec3(-556.0,0.0,0.0);
	Vec3 p0_3 = Vec3(556.0,0.0,559.2);
	recs.push_back(Rectangle(p0_3,a3,b3, Vec3(255,255,255)));
	//floor
	Vec3 b4 = Vec3(0.0,0.0,559.2);
	Vec3 a4 = Vec3(-556.0,0.0,0.0);
	Vec3 p0_4 = Vec3(556.0,0.0,0.0);
	recs.push_back(Rectangle(p0_4,a4,b4, Vec3(255,255,255)));
	//ceiling
	Vec3 a5 = Vec3(-556.0,0.0,0.0);
	Vec3 b5 = Vec3(0.0,0.0,-559.2);
	Vec3 p0_5 = Vec3(556.0,548.8,559.2);
	recs.push_back(Rectangle(p0_5,a5,b5, Vec3(255,255,255)));

	//short block top
	Vec3 a6 = Vec3(-50.0,0.0,160.0);
	Vec3 b6 = Vec3(160.0,0.0,50.0);
	Vec3 p0_6 = Vec3(130.0,165.0,65.0);
	recs.push_back(Rectangle(p0_6,a6,b6, Vec3(255,255,255)));

	//short block left 
	Vec3 b7 = Vec3(-50.0,0.0,160.0);
	Vec3 a7 = Vec3(0.0,165.0,0.0);
	Vec3 p0_7 = Vec3(290.0,0.0,115.0);
	recs.push_back(Rectangle(p0_7,a7,b7, Vec3(255,255,255)));

	//short block front
	Vec3 b8 = Vec3(160.0,0.0,50.0);
	Vec3 a8 = Vec3(0.0,165.0,0.0);
	Vec3 p0_8 = Vec3(130.0,0.0,65.0);
	recs.push_back(Rectangle(p0_8,a8,b8, Vec3(255,255,255)));

	//short block right
	Vec3 a9 = Vec3(-50.0,0.0,160.0);
	Vec3 b9 = Vec3(0.0,165.0,0.0);
	Vec3 p0_9 = Vec3(130.0,0.0,65.0);
	recs.push_back(Rectangle(p0_9,a9,b9, Vec3(255,255,255)));

	//short block back
	Vec3 a10 = Vec3(160.0,0.0,50.0);
	Vec3 b10 = Vec3(0.0,160.0,0.0);
	Vec3 p0_10 = Vec3(82.0,0.0,225.0);
	recs.push_back(Rectangle(p0_10,a10,b10, Vec3(255,255,255)));

	//tall block top
	Vec3 a11 = Vec3(-160.0,0.0,50.0);
	Vec3 b11 = Vec3(50.0,0.0,160.0);
	Vec3 p0_11 = Vec3(423.0,330.0,246.0);
	recs.push_back(Rectangle(p0_11,a11,b11, Vec3(255,255,255)));

	//tall block left
	Vec3 b12 = Vec3(50.0,0.0,160.0);
	Vec3 a12 = Vec3(0.0,330.0,0.0);
	Vec3 p0_12 = Vec3(423.0,0.0,246.0);
	recs.push_back(Rectangle(p0_12,a12,b12, Vec3(255,255,255)));

	//tall block back
	Vec3 b13 = Vec3(-160.0,0.0,50.0);
	Vec3 a13 = Vec3(0.0,330.0,0.0);
	Vec3 p0_13 = Vec3(473.0,0.0,406.0);
	recs.push_back(Rectangle(p0_13,a13,b13, Vec3(255,255,255)));

	//tall block right
	Vec3 b14 = Vec3(-50.0,0.0,-160.0);
	Vec3 a14 = Vec3(0.0,330.0,0.0);
	Vec3 p0_14 = Vec3(313.0,0.0,456.0);
	recs.push_back(Rectangle(p0_14,a14,b14, Vec3(255,255,255)));

	//tall block front
	Vec3 a15 = Vec3(-160.0,0.0,50.0);
	Vec3 b15 = Vec3(0.0,330.0,0.0);
	Vec3 p0_15 = Vec3(423.0,0.0,246.0);
	recs.push_back(Rectangle(p0_15,a15,b15, Vec3(255,255,255)));

	Viewport vp = Viewport(576,576,1);	
	ofstream img ("picture.ppm");

	img << "P3" << endl;
	img << vp.w << " " << vp.h << endl;
	img << "255" << endl;

	Vec3 w = dir * (-1);
	w = w / w.len();
	Vec3 u = up.cross(w);
	u = u / u.len();
	Vec3 v = w.cross(u);
	v = v / v.len();

	Ray ray = Ray();
	int n = sqrt(sample_per_pixel);
	multi_jitter(n);

	for (int y = vp.h-1; y >= 0; y--) {
		for (int x = 0; x < vp.w; x++) {
			//cout << "working on x:" << x << " y:" << y << endl;
			Vec3 pixel_value = Vec3(0,0,0);
			ray.o = eye;
			ray.d = vp.getPixelCenter(x,y,-eye.z);
			ray.d = u * ray.d.x + v * ray.d.y - w * ray.d.z;
			ray.d = ray.d / ray.d.len();

			vector<double> t_arr;
			for (int i = 0; i < recs.size(); i++) {
				double t;
				if (recs[i].hit(ray,t) == true) {
					t_arr.push_back(t);
				}
				else {
					t_arr.push_back(DBL_MAX);
				}
			}
			vector<double>::iterator it = min_element(t_arr.begin(),t_arr.end());
			int idx = it - t_arr.begin();
			double t = t_arr[idx];

			
			if (t != DBL_MAX) {
				if (idx == 0) {
					pixel_value = recs[0].material * ls;
					clamp(pixel_value);
					img << pixel_value <<endl;
				}
				else {
					Vec3 pt = ray.getPt(t);
					ray_count += area_light_num_rays;
					pixel_value = area_light_shade(pt,idx,area_light_num_rays);
					bool create = false;
					pixel_value = pixel_value + irradiance_caching(ray, pt, recs[idx].n,idx,ray_count,create,time);
					if (create) {
						clamp(pixel_value);
						img << pixel_value << endl;
						//img << "255 0 0" << endl;
					}
					else {
						clamp(pixel_value);
						img << pixel_value << endl;
					}
				}
			}
			/** only for irradiance caching **/
			else {
				img << "0 0 0" << endl;
			}
			/** Path tracing **/
		// 	Vec3 temp = Vec3(0,0,0);
		// 	for (int k = 0; k < sample_per_pixel; k++) {
		// 		ray_count ++;
		// 	// 	cout << ray_count << " " << x << " " << y << " " << k << endl;
				
		// 		ray.o = eye;
		// 		if (sample_per_pixel == 1)
		// 			ray.d = vp.getPixelCenter(x,y,-eye.z);
		// 		else
		// 			ray.d = vp.getPixelCenter(x + samples[2*k]-0.5,y+samples[2*k+1]-0.5,-eye.z);
		// 		ray.d = u * ray.d.x + v * ray.d.y - w * ray.d.z;
		// 		ray.d = ray.d / ray.d.len();
		// 	// 	//cout << samples[2*i] << endl;

		// 	// 	vector<double> t_arr;
		// 		for (int i = 0; i < recs.size(); i++) {
		// 			double t0;
		// 			if (recs[i].hit(ray,t0) == true) {
		// 				t_arr[i] = t0;
		// 			}
		// 			else {
		// 				t_arr[i] = DBL_MAX;
		// 			}
		// 		}
		// 		it = min_element(t_arr.begin(),t_arr.end());
		// 		idx = it - t_arr.begin();
		// 		t = t_arr[idx];

		// 		// if (t == DBL_MAX)
		// 		// 	continue;
		// 		if (t != DBL_MAX) {
		// 			Vec3 pt = ray.getPt(t);
		// 	// 		if (separate) {
		// 	// 			//cout << "not here" << endl;
		// 			// if (idx == 0) {
		// 			// 	pixel_value = pixel_value + recs[idx].material * ls;
		// 			// }
		// 			if (idx != 0) {
		// // 				//cout << indirect << endl;
		// 				temp = temp + path_shade(ray, 0, pt, idx, ray_count);
		// 			}
		// 		}
		// 	// 		else {
		// 	// 			if (idx == 0)
		// 	// 				pixel_value = pixel_value + recs[idx].material * ls;
		// 	// 			else {
		// 	// 				pixel_value = pixel_value + path_shade(ray, 0, pt, idx, ray_count);
		// 	// 			}
		// 	// 		}

		// 	// 		/** without shading **/
		// 	// 		//img << recs[idx].material << endl;
		// 	// 	}
		// 	// 	// if (pixel_value.x != 0 || pixel_value.y != 0 || pixel_value.z != 0)
		// 	// 	// cout << pixel_value << ", " << idx<< endl;
		// 	}
		// 	temp = temp / sample_per_pixel;
		// 	pixel_value = pixel_value + temp;
		// 	//cout << pixel_value << endl;
		// 	//cout << max(0,min(255, (int)(pixel_value.x))) << endl;
		// 	clamp(pixel_value);
		// 	img << pixel_value << endl;
		}
	}
	cout << (double)time / CLOCKS_PER_SEC << endl;
	//cout << (double)(clock() - startTime) / (CLOCKS_PER_SEC) << " seconds, " << ray_count << " rays" << endl;
	return 0;
}