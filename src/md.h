typedef struct {
	float x, y, z;
} float3;

typedef struct {
	int N;
	float3* r;
	float3* v;
	float3* f;
	float* m;
	int* type;
} MDSystem;
