#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include "DFT.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().

int main(int argc, char** argv)
{
    std::string op = std::string(argv[1]);

    Scene scene(784, 784);
    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);

    //add sphere
    //    Material* m = new Material(MICROFACET, Vector3f(0.0f));
    //    m->Ks = Vector3f(0.45f,0.45f,0.45f);
    //    m->Kd = Vector3f(0.3f,0.3f,0.3f);
    //    Sphere sphere1(Vector3f(150,100,300), 100, m);

    MeshTriangle floor("../models/cornellbox/floor.obj", white);
    MeshTriangle shortbox("../models/cornellbox/shortbox.obj", white);
    MeshTriangle tallbox("../models/cornellbox/tallbox.obj", white);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);

    scene.Add(&floor);
    //    scene.Add(&sphere1);
    scene.Add(&shortbox);
    scene.Add(&tallbox);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);

    scene.buildBVH();

    if(op.empty()||op[0]==48){
        Renderer r;

        auto start = std::chrono::system_clock::now();
        r.Render(scene);
        auto stop = std::chrono::system_clock::now();

        std::cout << "Render complete: \n";
        std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
        std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
        std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";
    }
    else{
        DFT dft;
        if(op[0]==49&&argc>2) dft.genFreImg(std::string(argv[2]));
        else if(op[0]==50&&argc>3) dft.genFreImg(std::string(argv[2]),std::string(argv[3]));
        else if(op[0]==51&&argc>2) dft.myFilter(std::string(argv[2]));
        else if(op[0]==52&&argc>2) dft.tilewiseSpec(std::string(argv[2]));
        else if(op[0]==53&&argc>2){
            Renderer r;
            auto start = std::chrono::system_clock::now();
            r.IMRender(scene, argv[2]);
            auto stop = std::chrono::system_clock::now();

            std::cout << "Render complete: \n";
            std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
            std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
            std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";
        }
    }


    return 0;
}