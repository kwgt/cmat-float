require 'matrix'

def print_c_source(m)
  print "    {\n"
  print "      #{m.row_size},\n"
  print "      #{m.row(0).size},\n"
  print "      (float[]) {\n"

  m.to_a.each {|row|
    tmp = row.inject([]) {|m, n| m << ("% 4d" % n)}
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
    matrix_info_t op1;
    float op2;
    matrix_info_t ans;
  } data[] = {
EOT

100.times {
  m = rand(15) + 1
  n = rand(15) + 1

  op1 = Matrix[*(Array.new(m) {Array.new(n) {rand(-10...+10)}})]
  op2 = rand(-10...+10)
  ans = op1 * op2

  print "  {\n"
  print_c_source(op1)
  print "    #{op2},\n"
  print_c_source(ans)
  print "  },\n"
}

print <<~EOT
 };
EOT
