#version 460
#extension GL_EXT_ray_tracing : enable

layout(location = 0) rayPayloadEXT bool payload;
layout(binding = 0, set = 0) uniform accelerationStructureEXT acc;
layout(binding = 1, set = 0, rgba8) uniform image2D image;
layout(binding = 2, set = 0) uniform Camera {
  vec3 corners[4];
  mat4 viewInverse;
} cam;

void main(void) {
  vec2  px        = vec2(gl_LaunchIDEXT.xy) + vec2(0.5);
  vec2  p         = px / vec2(gl_LaunchSizeEXT.xy);
  vec3  origin    = cam.viewInverse[3].xyz;
  vec3  target    = mix(mix(cam.corners[0], cam.corners[2], p.y), mix(cam.corners[1], cam.corners[3], p.y), p.x);
  vec4  direction = cam.viewInverse * vec4(normalize(target.xyz), 0.0);
  uint  rayFlags  = gl_RayFlagsOpaqueEXT | gl_RayFlagsSkipClosestHitShaderEXT;
  float tMin      = 0.1;
  float tMax      = 100.0;
  payload = true; // <- set payload to true as if we had hit the triangle (so, no need for a closest hit shader)
  traceRayEXT(
    acc,           // acceleration structure
    rayFlags,      // rayFlags
    0xFF,          // cullMask
    0,             // sbtRecordOffset // <- see comment [1] below
    0,             // sbtRecordStride // <- see comment [1] below
    0,             // missIndex
    origin,        // ray origin
    tMin,          // ray min range
    direction.xyz, // ray direction
    tMax,          // ray max range
    0              // payload (location = 0)
  );

  imageStore(image, ivec2(gl_LaunchIDEXT), 
    payload
      ? vec4(0.5, 0.6, 0.7, 1.0)
      : vec4(0.2, 0.3, 0.4, 1.0));
}