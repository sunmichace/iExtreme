<?php
function send_mes_exit($mes){
    header('Content-Type: application/json');
    echo json_encode($mes);
    exit;
}

function upload_file_size_is_zero($file_size){
    if($file_size==0){
        $mes = array(
            "success"=>0,
            "mes"=>"upload file is larger limitation",
            "data"=>""
        );
        send_mes_exit($mes);
    }
}

function valid_upload(){
    if(empty($_FILES)){
        $mes = array(
            "success"=>0,
            "mes"=>"upload file is empty",
            "data"=>""
        );
        send_mes_exit($mes);
    }
}

// function save_upload_file($des_dir, $user_hash){
//     move_uploaded_file($tmp_path, $save_file_path);
// };

function get_user_hash_flg(){
    $user_md5_str = $_SERVER['REMOTE_ADDR']."-".date('Y-m-d H:i:s').$_FILES["file"]["tmp_name"].$_FILES["file"]["size"];
    $user_md5_hash = md5($user_md5_str);
    return $user_md5_hash;
}

function select_mkdir($path){
    mkdir($path, 0777, false);
}

sleep(3);

#----------------------------
# valid
valid_upload();

$user_main_wd = getcwd()."/../model/imicrobes_user_pred_log";
$user_main_wd = str_replace('\\', '/', $user_main_wd);
$user_hash = get_user_hash_flg();
$user_dir = sprintf("%s/%s", $user_main_wd, $user_hash);
// file
// $pipe_file_path = sprintf("%s/%s/pipeline", $user_main_wd, $user_hash);
$raw_file_path = sprintf("%s/%s/raw", $user_main_wd, $user_hash);
$result_file_path = sprintf("%s/%s/result", $user_main_wd, $user_hash);
$user_cds_db = sprintf("%s/%s/cds_db", $user_main_wd, $user_hash);
select_mkdir($user_dir);
// select_mkdir($pipe_file_path);
select_mkdir($raw_file_path);
select_mkdir($result_file_path);
select_mkdir($user_cds_db);

// user upload file
$tmp_path = $_FILES["file"]["tmp_name"];
$file_name = "user_raw.fasta";
$cds_file_name = "cds_db.fna";
$file_size = $_FILES["file"]["size"];
upload_file_size_is_zero($file_size );

$user_raw_file_path = sprintf("%s/%s",$raw_file_path,$file_name);
$user_cds_db_file_path = sprintf("%s/%s",$user_cds_db,$cds_file_name);
// move
move_uploaded_file($tmp_path, $user_raw_file_path);
$cp_cds_db = sprintf("cp %s %s",$user_raw_file_path, $user_cds_db_file_path);
exec($cp_cds_db);

// select_mkdir($user_dir);
// save_upload_file($user_dir, $user_hash);

#----------------------------
# pre-processing
$input_path = $user_raw_file_path;
$outfile_path = $raw_file_path."/"."user_process.tsv";
$command = sprintf("python3 ../python/preprocess_data.py %s %s", $input_path, $outfile_path);
exec($command);

$mes = array(
    "success"=>1,
    "mes"=>"Upload file sucess, you can predict it now",
    "data"=>$user_hash
);
send_mes_exit($mes);
// print_r($_FILES);
?>