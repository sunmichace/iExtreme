



function waitBlock(){
    $(".wait-block").css({'transform': 'rotate(' + (++degree % 360) + 'deg)'});
    degree = degree>=360?0:degree
}

var degree = 0;
var order_id = null;
var is_uploaded = false;

$(function(){
    setInterval(waitBlock, 10);
})

function upload_file_event(event){
    event.preventDefault(); // 
}

// re-upload
function re_upload_event(){
    $("#file-input").on("click",function(){
        is_uploaded = false;
    });
}

function wait_block_tips(mes, is_open){
    if(is_open){
        $(".wait-block").show("slow");
        $("#wait-mes").show("slow");
        $("#wait-mes").text(mes);
    }else{
        $(".wait-block").hide("slow");
        $("#wait-mes").hide("slow");
        $("#wait-mes").text("");
    }
}

function showResult(data){
    $("#pred-results").show("slow");
    $("#vision-frame").show("slow");
    // $("#results-score").text(data["pred_score"].toFixed(5));
    // $("#VF-result").attr("href", data["vf_zip_url"]);
    $("#T-score").text(data["all_results"]["predict_score"][0]);
    $("#S-score").text(data["all_results"]["predict_score"][1]);
    $("#pH-score").text(data["all_results"]["predict_score"][2]);
    draw_TSpH(
        data["all_results"]["temp_x"],
        data["all_results"]["temp_y"],
        data["all_results"]["salinity_x"],
        data["all_results"]["salinity_y"],
        data["all_results"]["pH_x"],
        data["all_results"]["pH_y"],
    );
}

function tips_text_flash(mes){
    function hide_tips(){
        $("#tips-text").text("");
        $("#tips-text").hide("slow");
    }
    $("#tips-text").text(mes);
    $("#tips-text").show("slow");
    setTimeout(hide_tips,5000);
}

function ask_task(data){
    //
    if(data["status"]=="done"){
        is_uploaded = false;
        wait_block_tips("", false);
        showResult(data);
        tips_text_flash("Prediction complete");
        return
    }
    if(data["status"]=="error"){
        is_uploaded = false;
        wait_block_tips("", false);
        return
    }
    $.ajax({
        url:"/iExtreme_pred",
        type:"POST",
        data:{"req_type":"ask","order_id":order_id},
        dataType: 'json',
        success:function(res){
            // $(".elementor-counter-title").text(data["mes"]);
            wait_block_tips(res["mes"],true);
            setTimeout(function(){ask_task(res)},5000);
        },
        error:function(res){
            console.error(res);
        },
    });
}

function sumbit_ajax_init(){
    $("#btn-upload").click(function(){
        // $("#upload-results").toggle();
        // $("#file-upload").trigger('click');
        if(is_uploaded){
            // $("#tips-text").text("You have submitted");
            tips_text_flash("You have uploaded")
            return false;
        }
        $("#results-score").text("");
        $("#pred-results").hide();
        $("#vision-frame").hide();
        $("#VF-result").attr("href", "");
        is_uploaded =true;
        var formData = new FormData();
        var fileInput = $('#file-input')[0];
        var file = fileInput.files[0]; 
        formData.append('file', file);
        wait_block_tips("Uploading", true);
        $.ajax({
            url:"php/iExtreme_upload.php", // /web/php/ihalo_upload.php /ihalo_upload
            type:"post",
            data: formData,
            processData: false, 
            contentType: false,
            dataType:"json",
            success: function(res) {
                wait_block_tips("", false);
                tips_text_flash(res["mes"]);
                if(res["success"]==1){
                    $("#upload-results").text(res["mes"]);
                    $("#user-hash").attr("value", res["data"]);
                    $("#btn-pred").fadeIn("slow");
                }
            },
            error:function(res){
                console.error(res);
            }
        });
    });
    // predict btn
    $("#btn-pred").click(function(){
        $("#results-score").text("");
        $("#pred-results").hide();
        $("#btn-pred").hide("slow");
        var user_hash = $("#user-hash").attr("value");
        // console.log("user hash", user_hash);
        $.ajax({
            url:"/iExtreme_pred",
            type:"post",
            data: {"req_type":"submit", "user_hash":user_hash},
            success : function(data) {
                if(data["success"]==0){
                    is_uploaded = false;
                    $(".elementor-counter-title").text(data["mes"]);
                    $(".wait-block").css("display","none");
                    tips_text_flash(data["mes"]);
                    return
                }
                // success add in task line
                if(data["status"]=="add_success"){
                    order_id = data["order_id"];
                }
                if(data["status"]!="done"){
                    ask_task(data);
                }
            },
            error:function(res){
                console.error(res);
            }
        });
    });
}

$(function(){
    sumbit_ajax_init();
    re_upload_event();
})