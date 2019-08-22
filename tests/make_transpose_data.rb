require 'matrix'

def print_c_source(m)
  print "    {\n"
  print "      #{m.row_size},\n"
  print "      #{m.row(0).size},\n"
  print "      (float[]) {\n"

  m.to_a.each {|row|
    tmp = row.inject([]) {|m, n| m << ("% 3d" % n)}
    print "        #{tmp.join(",")},\n"
  }

  print "      },\n"
  print "    },\n"
end

print <<~EOT
  typedef struct {
    int rows;
    int cols;
    float* val;
  } matrix_info_t;

  static struct {
    matrix_info_t op;
    matrix_info_t ans;
  } data[] = {
EOT

100.times {
  r = rand(15) + 1
  c = rand(15) + 1

  op  = Matrix[*(Array.new(r) {Array.new(c) {rand(-10...+10)}})]
  ans = op.transpose

  print "  {\n"
  print_c_source(op)
  print_c_source(ans)
  print "  },\n"
}

print <<~EOT
 };
EOT
