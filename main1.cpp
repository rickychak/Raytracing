#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 2;

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& scene,
               int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.0f, 0.0f, 0.0f);
    Color r_color;

    Color ambient_color;
    glm::vec3 shadow_ray;
    glm::vec3 reflected_ray;
    float reflection;
    Material material = hit_record.material;

    ambient_color = material.k_a * material.ambient;
    color = ambient_color;
    //std::cout <<  hit_record.position[0] << hit_record.position[1] << std::endl;
    for (int i = 0; i < light_sources.size(); i++){
        shadow_ray = light_sources[i].position - hit_record.position;
        //std::cout <<  hit_record.position[0] << hit_record.position[1] << std::endl;
        
        if (glm::dot(hit_record.normal,shadow_ray) > 0){
            if (hittable_collection.Hit(shadow_ray,&hit_record)){
                //compute the second part of color
                reflection = pow(glm::dot(hit_record.reflection,hit_record.reflection),material.sh);
                color += light_sources[i].intensity * (material.k_d * material.diffuse + material.k_s * material.specular* reflection);      
            }
            else{

            }
        }
        
    }
    if (trace_depth < kMaxTraceDepth){
        if (material.k_s > 0){
            reflected_ray = hit_record.reflection;
            r_color= TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth++);
            color += material.k_s * r_color;
        }
    }

    if (color[0] >= 1){
        color[0] = 1;
    }
    if (color[1] >= 1){
        color[1] = 1;
    }
    if (color[2] >= 1){
        color[2] = 1;
    }

    return color;
}



Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record ;
    Color color(0.0f, 0.0f, 0.0f);

    
    if (hittable_collection.Hit(ray,&record)){

       //std::cout <<  record.position[0] << record.position[1] << std::endl;
       return Shade(light_sources,hittable_collection,record,trace_depth++);
    }
    return color;
}



int main() {
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("D://comp3271/Graphics_PA3_Release/Graphics_PA3_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/spheres.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++) {
                image[idx + i] = (uint8_t) (glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}

