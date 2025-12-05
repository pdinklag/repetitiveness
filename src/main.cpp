/**
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <sdsl/cst_sct3.hpp>

class Trie {
public:
    using Character = uint8_t;
    using NodeNumber = size_t;

private:
    static constexpr NodeNumber NIL = -1;
    static constexpr NodeNumber ROOT = 0;

    struct Node {
        Character label;
        NodeNumber first_child;
        NodeNumber next_sibling;
    } __attribute__((packed));

    std::vector<Node> nodes_;

    NodeNumber create_node(Character const label) {
        auto const x = nodes_.size();
        nodes_.push_back(Node{label, NIL, NIL});
        return x;
    }

public:
    Trie() {
        create_node(0); // root
    }

    bool try_get_child(NodeNumber const u, Character const c, NodeNumber& out) {
        NodeNumber prev = NIL;
        auto v = nodes_[u].first_child;
        while(v != NIL) {
            if(nodes_[v].label == c) {
                // move to front
                if(prev != NIL) {
                    nodes_[prev].next_sibling = nodes_[v].next_sibling;
                    nodes_[v].next_sibling = nodes_[u].first_child;
                    nodes_[u].first_child = v;
                }

                // return
                out = v;
                return true;
            }
            prev = v;
            v = nodes_[v].next_sibling;
        }
        return false;
    }

    NodeNumber insert_child(NodeNumber const u, Character const c) {
        auto const v = create_node(c);
        nodes_[v].next_sibling = nodes_[u].first_child;
        nodes_[u].first_child = v;
        return v;
    }

    NodeNumber root() const { return ROOT; }
};

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "usage: " << argv[0] << " <FILE> [prefix] [suffix-array] [lcp-array]" << std::endl;
        return -1;
    }

    size_t prefix = SIZE_MAX;
    if(argc >= 3) {
        prefix = std::atoll(argv[2]);
        if(prefix == 0) prefix = SIZE_MAX;
    }

    sdsl::cache_config cc;

    // load file
    std::cerr << "loading file ...";
    std::cerr.flush();

    sdsl::int_vector<8> text;
    {
        std::string s;
        {
            std::ifstream ifs(argv[1]);
            std::istreambuf_iterator<char> it(ifs);
            std::istreambuf_iterator<char> end;

            while(it != end && s.size() < prefix) {
                s.push_back(*it++);
            }
        }

        text = sdsl::int_vector<8>(s.size());
        for(size_t i = 0; i < s.size(); i++) {
            if(s[i] == 0 && i < s.size() - 1) {
                std::cerr << " failed -- the input file must not contain any zero bytes!" << std::endl;
                return -2;
            }
            text[i] = (unsigned char)s[i];
        }
        if(text[text.size()-1] != 0) sdsl::append_zero_symbol(text);
    }
    sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, cc);
    std::cerr << std::endl;

    auto const n = text.size();
    auto const actual_n = n - 1; // not taking into account the sentinel

    // construct SA
    sdsl::int_vector_buffer<> sa_buf;
    if(argc >= 4) {
        std::cerr << "loading SA ...";
        std::cerr.flush();

        sa_buf = sdsl::int_vector_buffer<>(argv[3]);
        cc.file_map[sdsl::conf::KEY_SA] = argv[3];
    } else {
        std::cerr << "computing SA ...";
        std::cerr.flush();

        sdsl::construct_sa<8>(cc);
        sa_buf = sdsl::int_vector_buffer<>(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
    }

    // cache SA in RAM
    sdsl::int_vector<> sa;
    sa.width(sa_buf.width());
    sa.resize(n);
    {
        for(size_t i = 0; i < n; i++) {
            sa[i] = sa_buf[i];
        }
    }

    std::cerr << std::endl;

    // output
    std::cout << "RESULT file=" << argv[1];

    // n
    std::cout << " n=" << actual_n; std::cout.flush();

    // alphabet and H0 entropy
    std::cout << " sigma="; std::cout.flush();
    size_t sigma = 0;
    double h0 = 0;
    {
        size_t hist[256];
        for(size_t c = 0; c < 256; c++) hist[c] = 0;

        for(size_t i = 0; i < actual_n; i++) ++hist[text[i]];

        for(size_t c = 0; c < 256; c++) {
            auto const nc = hist[c];
            if(nc) {
                ++sigma;
                h0 += (double(nc) / double(actual_n)) * std::log2(double(actual_n) / double(nc));
            }
        }
    }
    std::cout << sigma << " h0=" << h0; std::cout.flush();

    // BWT runs
    std::cout << " r="; std::cout.flush();
    size_t r = 0;
    {
        auto bwt = [&](size_t const i){
            auto const j = sa[i];
            return j > 0 ? text[j-1] : text[n-1];
        };
        
        uint8_t last = bwt(0);
        for(size_t i = 1; i < n; i++) {
            auto const c = bwt(i);
            if(last != 0 && c != last) ++r;
            last = c;
        }
    }
    std::cout << r; std::cout.flush();
    
    // LZ78
    std::cout << " z78="; std::cout.flush();
    size_t z78 = 0;
    {
        Trie trie;

        auto v = trie.root();
        for(size_t i = 0; i < actual_n; i++) {
            auto const c = text[i];
            if(!trie.try_get_child(v, c, v)) {
                trie.insert_child(v, c);
                v = trie.root();
                ++z78;
            }
        }
        if(v != trie.root()) ++z78; // final phrase
    }
    std::cout << z78;

    // LZ77
    std::cout << " z77="; std::cout.flush();
    size_t z77 = 0;
    {
        sdsl::construct_isa(cc);
        sdsl::int_vector_buffer<> isa(sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));
        
        for(size_t i = 0; i < actual_n;) {
            size_t const cur_pos = isa[i];

            // compute PSV and NSV
            auto lce = [&](size_t const i, size_t const j) {
                size_t l = 0;
                while(i + l < actual_n && j + l < actual_n && text[i + l] == text[j + l]) ++l;

                return l;
            };

            ssize_t psv_pos = (ssize_t)cur_pos - 1;
            while (psv_pos >= 0 && sa[psv_pos] > i) --psv_pos;
            size_t const psv_lcp = psv_pos >= 0 ? lce(i, sa[psv_pos]) : 0;

            size_t nsv_pos = cur_pos + 1;
            while(nsv_pos < actual_n && sa[nsv_pos] > i) ++nsv_pos;
            size_t const nsv_lcp = nsv_pos < actual_n ? lce(i, sa[nsv_pos]) : 0;

            // select maximum and advance
            auto const max_lcp = std::max(psv_lcp, nsv_lcp); // nb: may be zero
            auto const factor_len = std::max(size_t(1), max_lcp);
            i += factor_len;
            ++z77;
        }
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));
    }
    std::cout << z77; std::cout.flush();

    // delta -- courtesy of regindex/substring-complexity (MIT license)
    std::cout << " delta="; std::cout.flush();
    double delta = 0;
    {
        sdsl::int_vector_buffer<> lcp;
        if(argc >= 5) {
            lcp = sdsl::int_vector_buffer<>(argv[4]);
            cc.file_map[sdsl::conf::KEY_SA] = argv[4];
        } else {
            sdsl::construct_lcp_PHI<8>(cc);
            lcp = sdsl::int_vector_buffer<>(sdsl::cache_file_name(sdsl::conf::KEY_LCP, cc));
        }

        std::vector<uint32_t> dk(n, 0);
        for(size_t i = 1; i < n; i++) {
            dk[lcp[i]+1]++;
        }

        double x = dk[1];
        delta = x;
        for(size_t k = 2; k < n; k++) {
            x = x + dk[k] - 1;
            delta = std::max(delta, x / k);
        }
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_LCP, cc));
    }
    std::cout << std::fixed << delta; std::cout.flush();
    std::cout << std::endl;

    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, cc));

    return 0;
}
