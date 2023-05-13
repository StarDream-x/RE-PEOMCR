//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include "DFT.hpp"
#include <mutex>
#include <thread>


inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;

std::mutex mtx;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // change the spp value to change sample ammount
    int spp = 4;
    std::cout << "SPP: " << spp << "\n";

    //multi-thread
    int process = 0;
    const int thred = 20;
    int times = scene.height / thred;
    std::thread th[thred];

    auto castRayMultiThread = [&](uint32_t y_min, uint32_t y_max){
        for(uint32_t j = y_min; j < y_max;++j){
            int m = j * scene.width;
            for(uint32_t i = 0; i < scene.width; ++i){
                float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                          imageAspectRatio * scale;
                float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;
                Vector3f dir = normalize(Vector3f(-x, y, 1));
                for (int k = 0; k < spp; k++){
                    framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
                }
                m++;
            }
            mtx.lock();
            process++;
            UpdateProgress(process / (float)scene.height);
            mtx.unlock();
        }
    };

    for(int i = 0; i < thred; ++i)
        th[i] = std::thread(castRayMultiThread, i * times, (i + 1) * times);
    for(int i = 0; i< thred; ++i)
        th[i].join();
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}

double getLoss(Mat nImg, Mat GT){
    double sigma = sqrt(2.0 / acos(-1));
    Mat mImg,mGT;
    GaussianBlur(nImg, mImg, Size(3, 3), sigma, sigma);
//    GaussianBlur(GT,mGT,Size(3, 3), sigma, sigma);
//    nImg.copyTo(mImg);
    Mat difImg(nImg.size(),CV_32F);
    mImg.convertTo(mImg, CV_32F);
    GT.convertTo(GT, CV_32F);
    difImg = mImg - GT;
    return difImg.dot(difImg);
//    double res = mImg.at<float>(2,2) - GT.at<float>(2,2);
//    return res * res;
}

void Renderer::IMRender(const Scene &scene, std::string s) {
    int spp = 4, t = 10, zt = 15, md = 10;
    int znum[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    int er[] = {1,2,4,8};
    std::vector<Vector3f> framebuffer(scene.width * scene.height * spp);
    std::vector<Vector3f> nbuffer(scene.width * scene.height);
    GetVerticalSample(scene, spp, framebuffer);
//    printf("sample success!\n");

    Mat GT = imread(s);

    if(GT.empty()){
        std::cout<<"READ ERROR!"<<std::endl;
        exit(-1);
    }

//    printf("GOOD GT!\n");
    GT.convertTo(GT, CV_32FC3);
    Mat bgrGT[] = {Mat::zeros(GT.size(),CV_32F),
                    Mat::zeros(GT.size(),CV_32F),
                    Mat::zeros(GT.size(),CV_32F)};
    GT = GT / 255;
    split(GT, bgrGT);

    // 1. construct initial status
    for(int i = 0; i < scene.height * scene.width; ++i){
        nbuffer[i] = Vector3f(0,0,0);
        for(int k = 0; k < spp; ++k) if(er[k]&zt){
            nbuffer[i] += framebuffer[i * spp + k];
        }
        nbuffer[i] = nbuffer[i] / znum[zt];
        float temp = nbuffer[i].x;
        nbuffer[i].x = std::pow(clamp(0, 1, nbuffer[i].z), 0.6f);
        nbuffer[i].y = std::pow(clamp(0, 1, nbuffer[i].y), 0.6f);
        nbuffer[i].z = std::pow(clamp(0, 1, temp), 0.6f);
    }

    Mat imgMd(GT.size(),CV_32FC3,nbuffer.data());

    //2. split
    Mat bgrImg[] = {Mat::zeros(imgMd.size(),CV_32F),
                 Mat::zeros(imgMd.size(),CV_32F),
                 Mat::zeros(imgMd.size(),CV_32F)};
    split(imgMd,bgrImg);

    //3. IM
    int process = 0, mt = 0;
    for(int cn = 0; cn < 3; ++cn){
        for(int op = 0; op < t; ++op){
            for(int i = 0; i < scene.height; ++i){
                for(int j = 0; j < scene.width; ++j){
                    int nt = i * scene.width + j, lx = max(i-2-max(i + 3 - scene.height, 0),0), ly = max(j-2-max(j + 3 - scene.width, 0),0);
                    double lowLoss, gval;
                    for(int k = 1; k < 16; ++k){
                        float nval = 0.0f;
                        if(k&1) nval += framebuffer[nt*spp].at(2-cn);
                        if(k&2) nval += framebuffer[nt*spp+1].at(2-cn);
                        if(k&4) nval += framebuffer[nt*spp+2].at(2-cn);
                        if(k&8) nval += framebuffer[nt*spp+3].at(2-cn);
                        nval /= znum[k];
                        nval = std::pow(clamp(0,1,nval), 0.6f);
                        bgrImg[cn].at<float>(i,j) = nval;

//                        Mat m1(bgrImg[cn], Rect(lx, ly, 5, 5));
//                        Mat m1(bgrImg[cn], Range(lx,lx+5), Range(ly, ly+5));
//                        printf("md  (%d,%d)  (%d,%d)  ori=%.2f  ",i,j,lx+2,ly+2,bgrImg[cn].at<float>(lx+2,ly+2));
//                        printf("m1=%.2f   m2=%.2f\n",m1.at<float>(2,2),m2.at<float>(2,2));
                        double nLoss = getLoss(Mat(bgrImg[cn], Range(lx,lx+5), Range(ly,ly+5)),
                                               Mat(bgrGT[cn], Range(lx,lx+5), Range(ly,ly+5)));
//                        double nLoss = nval - bgrGT[cn].at<float>(i,j);
//                        nLoss = nLoss * nLoss;
                        if(k == 1 || lowLoss > nLoss){
                            lowLoss = nLoss;
                            gval = nval;
                            mt = 0;
                        }
                    }
                    bgrImg[cn].at<float>(i,j) = gval;
                    nbuffer[nt].set(cn,gval);
                }
                process++;
                UpdateProgress(process / (3.0f * scene.height * t));
            }
            ++mt;
            if(mt>md) break;
        }
    }

    //4. save
    FILE* fp = fopen("res.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        for(int j = 0; j < 3; ++j) color[j] = nbuffer[i].at(2-j) * 255;
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);
}

void Renderer::GetVerticalSample(const Scene& scene, int spp, std::vector<Vector3f> &framebuffer){

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // change the spp value to change sample ammount
    std::cout << "SPP: " << spp << "\n";

    //multi-thread
    int process = 0;
    const int thred = 20;
    int times = scene.height / thred;
    std::thread th[thred];

    auto castRayMultiThread = [&](uint32_t y_min, uint32_t y_max){
        for(uint32_t j = y_min; j < y_max;++j){
            int m = j * scene.width;
            for(uint32_t i = 0; i < scene.width; ++i){
                float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                          imageAspectRatio * scale;
                float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;
                Vector3f dir = normalize(Vector3f(-x, y, 1));
                for (int k = 0; k < spp; k++){
                    framebuffer[m * 4 + k] = scene.castRay(Ray(eye_pos, dir), 0);
                }
                m++;
            }
            mtx.lock();
            process++;
            UpdateProgress(process / (float)scene.height);
            mtx.unlock();
        }
    };

    for(int i = 0; i < thred; ++i)
        th[i] = std::thread(castRayMultiThread, i * times, (i + 1) * times);
    for(int i = 0; i< thred; ++i)
        th[i].join();
    UpdateProgress(1.f);
}