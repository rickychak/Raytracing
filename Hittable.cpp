#include "Hittable.h"
#include "math.h"

// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    float B = 2 * glm::dot(ray.o - o_, (ray.d));
    float C = glm::dot(ray.o-o_, ray.o-o_) - (pow(r_,2));
    float t;
    float determinant = pow(B, 2) - 4 * C;
    float D ;
    glm::vec3 tmp_Normal;
    if (determinant >= 0) {
        D = sqrt(determinant);
        float t1 = (-B + D) / 2;
        float t2 = (-B - D) / 2;
        if (t1 >= 0 && t2 < 0) {
            t = t1;
        }
        else if (t1 < 0 && t2 >= 0) {
            t = t2;
        }
        else if (t1 >= 0 && t2 >= 0) {
            if (t1 <= t2) {
                t = t1;
            }
            else {
                t = t2;
            }
        }
        else{
            return false;
        }
        if (t == 0) {
            return false;
        }
        tmp_Normal = glm::normalize(ray.At(t) - o_);
        hit_record->position = ray.At(t);
        hit_record->distance = glm::length(ray.d * t);
        hit_record->in_direction = ray.d;
        hit_record->material = material_;
        hit_record->normal = tmp_Normal;
        hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, tmp_Normal) * tmp_Normal);
        return true;
    }
    return false;
    
}

// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;
    glm::vec4 new_O = glm::vec4(ray.o, 1);
    glm::vec4 new_D = glm::vec4(ray.d, 0);
    float A = glm::dot(new_D, A_ * new_D);
    float B = 2*glm::dot(new_O, A_ * new_D);
    float C = glm::dot(new_O, A_ * new_O);
    float determinant = B * B - 4 * A * C;
    float t;
    float t1;
    float t2;
    if (determinant <0) {
        return false;
    }
    if (A != 0) {
        if (determinant == 0) {
            t = -B / (2 * A);
        }
        if (determinant > 0) {
            t1 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
            t2 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
            if (t1 > 0 && t2 > 0) {
                if (t1 >= t2) {
                    t = t2;
                }
                else {
                    t = t1;
                }
            }
            else {
                if (t1 >= 0) {
                    t = t1;
                }
                else {
                    t = t2;
                }
            }
            if (t == 0) {
                return false;
            }
        }
    }
    else {
        return false;
    }
    
    glm::vec3 tmp_Normal = glm::normalize((A_ + glm::transpose(A_)) * (new_O + new_D * t));
    glm::vec3 tmp_Reflection = ray.d - 2 * glm::dot(ray.d, tmp_Normal) * tmp_Normal;
    hit_record->position = ray.At(t);
    hit_record->distance = sqrt(dot(t*ray.d,t*ray.d));
    hit_record->in_direction = ray.d;
    hit_record->material = material_;
    hit_record->normal = tmp_Normal;
    hit_record->reflection =glm::normalize( tmp_Reflection);
    return true;

}


// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;
    glm::vec3 tmp_Plane = -glm::cross((c_ - a_), (b_ - a_));
    glm::vec3 unit_plane_normal = glm::normalize(tmp_Plane);
    float total_Area = glm::length(tmp_Plane) / 2;
    float t = (dot(a_, unit_plane_normal) - dot(ray.o, unit_plane_normal))/dot(unit_plane_normal, ray.d);
    if (t <= 0) {
        return false;
    }
    glm::vec3 intersect_Point = ray.At(t);
    
    float area_1 = glm::length(glm::cross(a_ - intersect_Point, b_ - intersect_Point))/2;
    float area_2 = glm::length(glm::cross((b_-intersect_Point ), (c_-intersect_Point ))) / 2;
    float area_3 = glm::length(glm::cross((c_-intersect_Point ), (a_-intersect_Point ))) / 2;
    

    
    glm::vec3 cross1 = glm::cross((b_ - a_), (intersect_Point - a_));
    glm::vec3 cross2 = glm::cross((intersect_Point - a_), (c_ - a_));
    glm::vec3 cross3 = glm::cross((intersect_Point - c_), (b_ - c_));

    
   
    float alpha = area_1 / total_Area;
    float beta = area_2 / total_Area;
    float gamma = 1-alpha-beta;

    
    if (!(glm::dot(cross1, unit_plane_normal) >= 0 && glm::dot(cross2, unit_plane_normal) >= 0 && glm::dot(cross3, unit_plane_normal) >= 0)) {
        return false;
    }
    else{
        hit_record->position = ray.At(t);
        hit_record->distance = glm::length(t*ray.d);
        hit_record->in_direction = ray.d;
        
        if (phong_interpolation_) {
            
            glm::vec3 phong_Normal = (beta * n_a_ + gamma * n_b_ +alpha * n_c_);
            glm::vec3 unit_Normal = glm::normalize(phong_Normal);
            hit_record->normal = unit_Normal;
            glm::vec3 tmp_Reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, unit_Normal) * unit_Normal);
            hit_record->reflection = tmp_Reflection;
        }
        else {
            hit_record->normal = unit_plane_normal;
            glm::vec3 tmp_Reflection = ray.d - 2 * glm::dot(ray.d, unit_plane_normal) * unit_plane_normal;
            hit_record->reflection = glm::normalize(tmp_Reflection);

        }
        // no need to set material in this function
        return true;
    }
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret) {
        hit_record->material = material_;
    }
    return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
           const Material& material,
           bool phong_interpolation):
           ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto& face : f_ind_) {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto& vertex_normal : vertex_normals_) {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto& face : f_ind_) {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min( 1e5f,  1e5f,  1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto& vertex : vertices_) {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++) {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray& ray, HitRecord *hit_record) const {
    const bool brute_force = false;
    if (brute_force) {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_) {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record)) {
                if (curr_hit_record.distance < min_dist) {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f) {
            hit_record->material = material_;
            return true;
        }
        return false;
    } else {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret) {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
    for (size_t idx : face) {
        const auto& pt = vertices_[idx];
        for (int i = 0; i < 3; i++) {
            if (pt[i] < bbox_min[i] + 1e-6f) return false;
            if (pt[i] > bbox_max[i] - 1e-6f) return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
    float t_min = -1e5f;
    float t_max =  1e5f;

    for (int i = 0; i < 3; i++) {
        if (glm::abs(ray.d[i]) < 1e-6f) {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
                t_min =  1e5f;
                t_max = -1e5f;
            }
        }
        else {
            if (ray.d[i] > 0.f) {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
    const Point& bbox_min = u->bbox_min;
    const Point& bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++) {
        for (size_t b = 0; b < 2; b++) {
            for (size_t c = 0; c < 2; c++) {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
                    if (u->childs[child_idx] == nullptr) {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode* child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs) {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
        return false;
    }
    float distance = 1e5f;
    for (const auto& face_idx : u->face_index) {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto& child : u->childs) {
        if (child == nullptr) {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord *hit_record) const {
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_) {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < min_dist) {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}