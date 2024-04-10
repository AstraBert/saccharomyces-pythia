import gradio as gr
import os
import time
from transformers import pipeline
from predict import *
from load_model import *

def print_like_dislike(x: gr.LikeData):
    print(x.index, x.value, x.liked)

def add_message(history, message):
    if len(message["files"]) > 0:
        history.append((message["files"], None))
    if message["text"] is not None and message["text"] != "":
        history.append((message["text"], None))
    return history, gr.MultimodalTextbox(value=None, interactive=False)


def bot(history):
    global tsk
    if type(history[-1][0]) != tuple:
        try:
            pipe = pipeline("text-generation", tokenizer=tokenizer, model=model)
            response = pipe(history[-1][0])[0]
            response = response["generated_text"]
            history[-1][1] = ""
            for character in response:
                history[-1][1] += character
                time.sleep(0.05)
                yield history
        except Exception as e:
            response = f"Sorry, the error '{e}' occured while generating the response; check [troubleshooting documentation](https://astrabert.github.io/everything-rag/#troubleshooting) for more"
    if type(history[-1][0]) == tuple:
        filelist = []
        for i in history[-1][0]:
            filelist.append(i)
        if len(filelist) > 1:
            finalfasta = merge_fastas(filelist)
        else:
            finalfasta = filelist[0]
        response = predict_genes(finalfasta)
        history[-1][1] = ""
        for character in response:
            history[-1][1] += character
            time.sleep(0.05)
            yield history

with gr.Blocks() as demo:
    chatbot = gr.Chatbot(
        [[None, " Welcome to Saccharomyces-Pythia, your helpful assistant for all things Saccharomyces cerevisiae! I am here to provide you with fascinating facts about this important model organism, as well as aid in the prediction of open reading frames (ORFs) and their corresponding types from any S. cerevisiae genetic sequence you may have. Simply upload your FASTA file, and let me work my magic. Rest assured, accuracy and efficiency are at the core of my design. Prepare to be enlightened on the wonders of yeast genomics and beyond. Let's get started!"]],
        label="Saccharomyces-Pythia",
        elem_id="chatbot",
        bubble_full_width=False,
    )

    chat_input = gr.MultimodalTextbox(interactive=True, file_types=["pdf"], placeholder="Enter message or upload file...", show_label=False)

    chat_msg = chat_input.submit(add_message, [chatbot, chat_input], [chatbot, chat_input])
    bot_msg = chat_msg.then(bot, chatbot, chatbot, api_name="bot_response")
    bot_msg.then(lambda: gr.MultimodalTextbox(interactive=True), None, [chat_input])

    chatbot.like(print_like_dislike, None, None)
    clear = gr.ClearButton(chatbot)

demo.queue()
if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", share=False)

	