
#[allow(dead_code)]
const fn assert_tuple_is_packed() {
    assert!(std::mem::size_of::<(f64, f64)>() == 16);
}
const _: () = assert_tuple_is_packed();

pub extern "C" fn iw_pip_inside_simd(vertices: *const f64, count: usize, px: f64, py: f64) -> bool {
    let p = (px, py);
    unsafe {
        // transmute is safe because assert_tuple_is_packed is valid.
        let pointpointer = std::mem::transmute::<_, *const (f64, f64)>(vertices);
        let vertex_slice: &[(f64, f64)] = std::slice::from_raw_parts(pointpointer, count);
        crate::inside_simd(vertex_slice, &p)
    }
}


pub extern "C" fn iw_pip_edge_tree_create(vertices: *const f64, count: usize) -> *const crate::edge_tree::EdgeTree {
    unsafe {
        // transmute is safe because assert_tuple_is_packed is valid.
        let pointpointer = std::mem::transmute::<_, *const (f64, f64)>(vertices);
        let vertex_slice: &[(f64, f64)] = std::slice::from_raw_parts(pointpointer, count);
        let b = Box::new(crate::edge_tree::EdgeTree::new(vertex_slice));
        Box::leak(b)
    }
}

pub extern "C" fn iw_pip_edge_tree_test(tree: *const crate::edge_tree::EdgeTree, px: f64, py: f64) ->  bool{
    unsafe {
        (*tree).inside(&(px, py))
    }
}

pub extern "C" fn iw_pip_edge_tree_free(tree: *mut crate::edge_tree::EdgeTree) {
    unsafe {
        let boxed_thing = Box::from_raw(tree);
        drop(boxed_thing);
    }
}


