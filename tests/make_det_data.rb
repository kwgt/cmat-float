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
    float ans;
  } data[] = {
EOT

100.times { |i|
  n = rand(15) + 1

  op  = Matrix[*(Array.new(n) {Array.new(n) {rand(-10...+10)}})]
  ans = op.det

  print "  // #{i}\n"
  print "  {\n"
  print_c_source(op)
  print "    #{ans}\n"
  print "  },\n"
}

print <<~EOT
 };
EOT
