require 'matrix'
require 'prime'

class Integer
  def divisor_list
    return [1] if self == 1

    # 数を素因数分解
    pf = Prime.prime_division(self)
      
    # 素因数毎に因数要素のリストに展開
    pf.map! {|e| Array.new(e[1]+1) {|i| e[0]**i}}

    # 因数要素の組み合わせを作成
    fl = pf.inject{|m, n| m.product(n)}

    # 因数要素を掛け合わせて、因数のリストに変換
    fl.map! {|e| [e].flatten.inject(&:*)}

    return fl.sort
  end
end

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
    matrix_info_t op1;
    matrix_info_t op2;
    float ans;
  } data[] = {
EOT

100.times {
  begin
    s = rand(60) + 1
  end  while Prime.prime?(s)

  fl = s.divisor_list
  n  = fl[rand(fl.size)]
  m  = fl[rand(fl.size)]

  op1 = Matrix[*(Array.new(n) {Array.new(s / n) {rand(-10...+10)}})]
  op2 = Matrix[*(Array.new(m) {Array.new(s / m) {rand(-10...+10)}})]

  v1  = Vector[*op1.to_a.flatten]
  v2  = Vector[*op2.to_a.flatten]

  print "  {\n"
  print_c_source(op1)
  print_c_source(op2)
  print "    #{v1.inner_product(v2)}\n"
  print "  },\n"
}

print <<~EOT
 };
EOT
