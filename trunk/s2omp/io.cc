#include <iostream>
#include <fstream>

#include "io.h"

#include "pixel.h"
#include "field_pixel-inl.h"
#include "tree_pixel.h"
#include "point.h"
#include "pixel_union.h"
#include "field_union.h"
#include "tree_union.h"

namespace s2omp {

bool io::write_ascii(const point_vector& points, const string& file_name) {
  std::fstream fs(file_name, std::fstream::out);
  if (!fs) {
    return false;
  }

  for (int k = 0; k < points.size(); k++) {
    fs << points[k].to_pixel() << " " << points[k].weight() << "\n";
  }
  fs.close();

  return true;
}

bool io::write_ascii(const pixel_vector& pixels, const string& file_name) {
  std::fstream fs(file_name, std::fstream::out);
  if (!fs) {
    return false;
  }

  for (int k = 0; k < pixels.size(); k++) {
    fs << pixels[k] << "\n";
  }
  fs.close();

  return true;
}

bool io::write_ascii(const pixel_union& pix_union, const string& file_name) {
  std::fstream fs(file_name, std::fstream::out);
  if (!fs) {
    return false;
  }

  for (pixel_iterator iter = pix_union.begin(); iter != pix_union.end(); ++iter) {
    fs << *iter << "\n";
  }
  fs.close();

  return true;
}

bool io::write_pb(const point_vector& points, const string& output_file) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  PointVectorProto point_vector_proto;
  for (int k = 0; k < points.size(); k++) {
    PointProto* point = point_vector_proto.add_point();
    point->set_id(points[k].id());
    point->set_weight(points[k].weight());
  }

  std::fstream fs(output_file, std::ios::out | std::ios::trunc
      | std::ios::binary);
  return point_vector_proto.SerializeToOstream(&fs);
}

bool io::write_pb(const pixel_vector& pixels, const string& output_file) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  PixelVectorProto pixel_vector_proto;
  for (int k = 0; k < pixels.size(); k++) {
    PixelProto* pixel = pixel_vector_proto.add_pixel();
    pixel->set_id(pixels[k].id());
  }

  std::fstream fs(output_file, std::ios::out | std::ios::trunc
      | std::ios::binary);
  return pixel_vector_proto.SerializeToOstream(&fs);
}

bool io::write_pb(const pixel_union& pix_union, const string output_file) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  PixelVectorProto pixel_vector_proto;
  for (pixel_iterator iter = pix_union.begin(); iter != pix_union.end(); ++iter) {
    PixelProto* pixel = pixel_vector_proto.add_pixel();
    pixel->set_id(iter->id());
  }

  std::fstream fs(output_file, std::ios::out | std::ios::trunc
      | std::ios::binary);
  return pixel_vector_proto.SerializeToOstream(&fs);
}

bool io::write_pb(const field_union& field, const string& output_file) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  FieldUnionProto field_union_proto;
  switch (field.type()) {
  case field_union::SCALAR_FIELD:
    field_union_proto.set_field_type(FieldUnionProto::SCALAR_FIELD);
    break;
  case field_union::DENSITY_FIELD:
    field_union_proto.set_field_type(FieldUnionProto::DENSITY_FIELD);
    break;
  case field_union::SAMPLED_FIELD:
    field_union_proto.set_field_type(FieldUnionProto::SAMPLED_FIELD);
    break;
  }

  for (field_const_iterator iter = field.begin(); iter != field.end(); ++iter) {
    PixelProto* pixel = field_union_proto.add_pixel();
    pixel->set_id(iter->id());
    pixel->set_weight(iter->weight());
    pixel->set_intensity(iter->intensity());
    pixel->set_n_points(iter->n_points());
  }

  std::fstream fs(output_file, std::ios::out | std::ios::trunc
      | std::ios::binary);
  return field_union_proto.SerializeToOstream(&fs);
}

bool io::write_pb(const tree_union& t, const string& output_file) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  TreeUnionProto tree_union_proto;
  tree_union_proto.set_level(t.level());
  tree_union_proto.set_node_capacity(t.node_capacity());

  point_vector points;
  t.copy_points(&points);
  for (int k = 0; k < points.size(); k++) {
    PointProto* point = tree_union_proto.add_point();
    point->set_id(points[k].id());
    point->set_weight(points[k].weight());
  }

  std::fstream fs(output_file, std::ios::out | std::ios::trunc
      | std::ios::binary);
  return tree_union_proto.SerializeToOstream(&fs);
}

bool io::read_ascii(const string& input_file, point_vector* points) {
  if (!points->empty())
    points->clear();

  std::fstream fs(input_file, std::ios::in);
  if (!fs) {
    return false;
  }

  while (!fs.eof()) {
    pixel pix;
    double weight;
    fs >> pix >> weight;

    if (!fs.eof()) {
      points->push_back(point(pix.id(), weight));
    }
  }

  return true;
}

bool io::read_ascii(const string& input_file, pixel_vector* pixels) {
  if (!pixels->empty())
    pixels->clear();

  std::fstream fs(input_file, std::ios::in);
  if (!fs) {
    return false;
  }

  while (!fs.eof()) {
    pixel pix;
    fs >> pix;

    if (!fs.eof()) {
      pixels->push_back(pix);
    }
  }

  return true;
}

bool io::read_ascii(const string& input_file, pixel_union* pix_union) {
  pixel_vector pixels;

  if (!read_ascii(input_file, &pixels)) {
    return false;
  }

  pix_union->init(pixels);

  return true;
}

bool io::read_pb(const string& input_file, point_vector* points) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (!points->empty())
    points->clear();

  PointVectorProto proto;
  std::fstream fs(input_file, std::ios::in | std::ios::binary);
  if (!fs || !proto.ParseFromIstream(&fs)) {
    return false;
  }

  for (int k = 0; k < proto.point_size(); k++) {
    points->push_back(to_point(proto.point(k)));
  }

  return true;
}

bool io::read_pb(const string& input_file, pixel_vector* pixels) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (!pixels->empty())
    pixels->clear();

  PixelVectorProto proto;
  std::fstream fs(input_file, std::ios::in | std::ios::binary);
  if (!fs || !proto.ParseFromIstream(&fs)) {
    return false;
  }

  for (int k = 0; k < proto.pixel_size(); k++) {
    pixels->push_back(to_pixel(proto.pixel(k)));
  }

  return true;
}

bool io::read_pb(const string& input_file, pixel_union* pix_union) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  pixel_vector pixels;
  if (!read_pb(input_file, &pixels)) {
    return false;
  }

  pix_union->init(pixels);

  return true;
}

bool io::read_pb(const string& input_file, field_union* s) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  FieldUnionProto proto;
  std::fstream fs(input_file, std::ios::in | std::ios::binary);
  if (!fs || !proto.ParseFromIstream(&fs)) {
    return false;
  }

  field_union::FieldType type;
  switch (proto.field_type()) {
  case FieldUnionProto::SCALAR_FIELD:
    type = field_union::SCALAR_FIELD;
    break;
  case FieldUnionProto::DENSITY_FIELD:
    type = field_union::DENSITY_FIELD;
    break;
  case FieldUnionProto::SAMPLED_FIELD:
    type = field_union::SAMPLED_FIELD;
    break;
  }

  field_vector pixels;
  for (int k = 0; k < proto.pixel_size(); k++) {
    pixels.push_back(to_field_pixel(proto.pixel(k)));
  }

  return s->init(pixels, type);
}

bool io::read_pb(const string& input_file, tree_union* t) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  TreeUnionProto proto;
  std::fstream fs(input_file, std::ios::in | std::ios::binary);
  if (!fs || !proto.ParseFromIstream(&fs)) {
    return false;
  }

  t->init(proto.level(), proto.node_capacity());

  for (int k = 0; k < proto.point_size(); k++) {
    if (!t->add_point(to_point(proto.point(k)))) {
      return false;
    }
  }

  return true;
}

point io::to_point(const PointProto& proto) {
  return point(proto.id(), proto.weight());
}

pixel io::to_pixel(const PixelProto& proto) {
  return pixel(proto.id());
}

field_pixel io::to_field_pixel(const PixelProto& proto) {
  return field_pixel(proto.id(), proto.intensity(),
      proto.weight(), proto.n_points());
}

} // end namespace s2omp
