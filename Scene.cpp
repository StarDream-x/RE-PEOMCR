//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here

    Intersection intersection = intersect(ray);
    if(!intersection.happened) return Vector3f(0.0f, 0.0f, 0.0f);
    if(intersection.m->hasEmission()) return intersection.m->getEmission();

    Vector3f L_dir(0.0f,0.0f,0.0f);
    Vector3f L_indir(0.0f,0.0f,0.0f);
    Intersection inter;
    float pdf_light;
    sampleLight(inter, pdf_light);
    Vector3f dir = inter.coords - intersection.coords;
    float dis = dotProduct(dir,dir);
    Intersection ckInter = intersect(Ray(intersection.coords,dir.normalized()));
    if(ckInter.distance - dir.norm() > -0.005f)
    {
        dir = dir.normalized();
        L_dir = inter.emit * intersection.m->eval(ray.direction,dir,intersection.normal)
                * dotProduct(dir,intersection.normal) * dotProduct(-dir, inter.normal)
                / dis / pdf_light;
    }

    float rd = get_random_float();
    if(rd < RussianRoulette){
        Vector3f wi = intersection.m->sample(ray.direction,intersection.normal);
        wi = wi.normalized();
        inter = intersect(Ray(intersection.coords,wi));
        if((inter.happened)&&(!inter.m->hasEmission())){
            L_indir = castRay(Ray(intersection.coords,wi),depth+1)
                    * intersection.m->eval(ray.direction, wi, intersection.normal)
                    * dotProduct(wi, intersection.normal)
                    / intersection.m->pdf(ray.direction,wi,intersection.normal)
                    / RussianRoulette;
        }
    }

    return L_dir + L_indir;
}