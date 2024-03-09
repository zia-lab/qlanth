(*-----------------------------------------------------------------------+
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~+----------------------------------------------------------------+~~~|
|~~~|                     __           __   _                        |~~~|
|~~~|               _____/ /___ ______/ /__(_)___  ____ _            |~~~|
|~~~|              / ___/ / __ `/ ___/ //_/ / __ \/ __ `/            |~~~|
|~~~|             (__  ) / /_/ / /__/ ,< / / / / / /_/ /             |~~~|
|~~~|            /____/_/\__,_/\___/_/|_/_/_/ /_/\__, /              |~~~|
|~~~|                                           /____/               |~~~|
|~~~|                                                                |~~~|
|~~~+----------------------------------------------------------------+~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~+--------------------------------------------------------------+~~~~|
|~~~~|                                                              |~~~~|
|~~~~|   This package is useful for posting messages, files, and    |~~~~|
|~~~~|                    expressions to Slack.                     |~~~~|
|~~~~|   Good for monitoring processes that take a long time to     |~~~~|
|~~~~|              run and for keeping lab records.                |~~~~|
|~~~~|                                                              |~~~~|
|~~~~+--------------------------------------------------------------+~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
+-----------------------------------------------------------------------*)

keys          = Import["keys.m"]
slackToken    = keys["slackToken"];
slackUserName = keys["slackUserName"];

PostMessageToSlack::usage="PostMessageToSlack[text, slackChannel, Options] posts a message to a Slack channel. Options include \"threadTS\" -> None, which specifies the thread to post the message to, if None, the message is posted as the head of a new thread.";
Options[PostMessageToSlack] = {"threadTS" -> None};
PostMessageToSlack[text_, slackChannel_, OptionsPattern[]] :=
    Module[{params, response},
        params = {
            "token" -> slackToken,
            "channel" -> slackChannel,
            "text" -> text,
            "username" -> slackUserName
        };
        threadTS = OptionValue["threadTS"];
        If[threadTS =!= None, AppendTo[params, "thread_ts" -> threadTS]];

        response = URLExecute[HTTPRequest["https://slack.com/api/chat.postMessage", 
            <|
                "Method" -> "POST",
                "Body" -> params
            |>],
            "RawJSON"
        ];

        If[!KeyExistsQ[response, "ok"] || ! response["ok"],
            Print["Error in request: ", response];
            Return[$Failed]
        ];
        Return[response];
    ];

PostFileToSlack::usage="PostFileToSlack[text, filePath, slackChannel, Options] posts a file to a Slack channel. Options include \"fileType\" -> None, \"title\" -> None, \"threadTS\" -> None, which specifies the thread to post the message to, if None, the message is posted as the head of a new thread.";
Options[PostFileToSlack] = {"fileType" -> None, "title" -> None, "threadTS" -> None};
PostFileToSlack[text_, filePath_, slackChannel_, OptionsPattern[]] := 
    Module[{params, response, fileName, fileBytes},  
        threadTS = OptionValue["threadTS"];
        fileType = OptionValue["fileType"];
        title = OptionValue["title"];
        fileName = FileNameTake[filePath];
        params = {
            "token" -> slackToken,
            "channels" -> slackChannel, 
            "filename" -> fileName,
            "initial_comment" -> text, 
            "file" -> File[fileName]};
        If[fileType =!= None, AppendTo[params, "filetype" -> fileType]];
        If[title =!= None,    AppendTo[params, "title" -> title]];
        If[threadTS =!= None, AppendTo[params, "thread_ts" -> threadTS]];
        response = URLExecute[HTTPRequest["https://slack.com/api/files.upload", 
            <|
                "Method" -> "POST", 
                "Body" -> params
            |>
        ],
        "RawJSON"
        ];
        If[! KeyExistsQ[response, "ok"] || ! response["ok"], 
            Print["Error in request: ", response];
            Return[$Failed]
        ];
        Return[response];
    ];

PostPdfToSlack::usage="PostPDFToSlack[text, expression, slackChannel, Options] posts expression to the given Slack channel acccompanied by the given text. The expression is exported to an intermediate pdf and that is the file posted. Options include  \"title\" -> None, this specifies the title of the file,
\"threadTS\" is another option which specifies the thread to post the message to, if None, the message is posted as the head of a new thread.";
Options[PostPdfToSlack] = {"title" -> None, "threadTS" -> None};
PostPdfToSlack[text_, graph_, slackChannel_, OptionsPattern[]] :=  
    Module[{params, response, fileName, fileBytes}, 
        threadTS = OptionValue["threadTS"];
        title    = OptionValue["title"];
        fileName = FileNameTake[filePath];
        filePath = "temp_graph.pdf";
        Export[filePath, graph];
        params = {
            "token" -> slackToken,
            "channels" -> slackChannel, 
            "filename" -> fileName,
            "initial_comment" -> text, 
            "file" -> File[fileName]};
        If[title    =!= None, AppendTo[params, "title" -> title]];
        If[threadTS =!= None, AppendTo[params, "thread_ts" -> threadTS]];
        response = URLExecute[HTTPRequest["https://slack.com/api/files.upload",
            <|
                "Method" -> "POST", 
                "Body" -> params
            |>],
            "RawJSON"
        ];
        If[! KeyExistsQ[response, "ok"] || ! response["ok"], 
            Print["Error in request: ", response];
            Return[$Failed]
        ];
        response
    ];
